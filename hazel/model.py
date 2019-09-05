from hazel.chromosphere import Hazel_atmosphere
from hazel.photosphere import SIR_atmosphere
from hazel.parametric import Parametric_atmosphere
from hazel.stray import Straylight_atmosphere
from hazel.configuration import Configuration
from hazel.io import Generic_output_file
from collections import OrderedDict
from hazel.codes import hazel_code, sir_code
from hazel.spectrum import Spectrum
from hazel.transforms import transformed_to_physical, physical_to_transformed, jacobian_transformation
import hazel.util
import numpy as np
import copy
import os
from pathlib import Path
import scipy.stats
import scipy.special
import scipy.signal
import scipy.linalg
import scipy.optimize
import warnings
import logging

# from ipdb import set_trace as stop

__all__ = ['Model']

class Model(object):
    def __init__(self, config=None, working_mode='synthesis', verbose=0, debug=False, rank=0, randomization=None):

        np.random.seed(123)

        if (rank != 0):
            return
        
        self.photospheres = []
        self.chromospheres = []
        self.chromospheres_order = []
        self.atmospheres = {}
        self.order_atmospheres = []        
        self.straylight = []
        self.parametric = []
        self.spectrum = []
        self.configuration = None
        self.n_cycles = 1
        self.spectrum = {}
        self.topologies = []
        self.straylights = []
        self.working_mode = working_mode
        self.pixel = 0
        self.debug = debug

        self.epsilon = 1e-2
        self.step_limiter_inversion = 1.0
        self.backtracking = 'brent'
        
        self.verbose = verbose
        
        self.logger = logging.getLogger("model")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []

        ch = logging.StreamHandler()
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        # Set randomization
        if (randomization is None):
            self.n_randomization = 1        
        else:
            self.n_randomization = randomization
        
        if (config is not None):
            if (self.verbose >= 1):
                self.logger.info('Using configuration from file : {0}'.format(config))
            self.configuration = Configuration(config)

            self.use_configuration(self.configuration.config_dict)
        
        # Initialize pyhazel
        hazel_code._init()        

    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d:
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d:
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)
            
    def use_configuration(self, config_dict):
        """
        Use a configuration file

        Parameters
        ----------
        config_dict : dict
            Dictionary containing all the options from the configuration file previously read

        Returns
        -------
        None
        """        

        # Deal with the spectral regions        
        tmp = config_dict['spectral regions']


        # Output file
        self.output_file = config_dict['working mode']['output file']

        # Backtracking mode
        if ('backtracking' in config_dict['working mode']):
            self.backtracking = config_dict['working mode']['backtracking']
        else:
            self.backtracking = 'brent'
        if (self.verbose >= 1):
            self.logger.info('Backtracking mode : {0}'.format(self.backtracking))
        
        
        # Working mode
        # self.working_mode = config_dict['working mode']['action']

        # Add spectral regions
        for key, value in config_dict['spectral regions'].items():
            self.add_spectral(value)

        # Set number of cycles if present
        if (self.working_mode == 'inversion'):
            if ('number of cycles' in config_dict['working mode']):
                if (config_dict['working mode']['number of cycles'] != 'None'):
                    self.n_cycles = int(config_dict['working mode']['number of cycles'])
                    if (self.verbose >= 1):
                        self.logger.info('Using {0} cycles'.format(self.n_cycles))

        # Set number of maximum iterations
        if ('maximum iterations' in config_dict['working mode']):
            if (config_dict['working mode']['number of cycles'] != 'None'):
                self.max_iterations = int(config_dict['working mode']['maximum iterations'])
                if (self.verbose >= 1):
                    self.logger.info('Using {0} max. iterations'.format(self.max_iterations))
            else:
                self.max_iterations = 10
        else:
            self.max_iterations = 10

        # Randomization
        if (self.verbose >= 1):
            if (self.n_randomization == 1):
                self.logger.info('Not using randomizations')
            else:
                self.logger.info('Using a maximum of {0} randomizations'.format(self.n_randomization))

        # Set number of maximum iterations
        if ('relative error' in config_dict['working mode']):
            if (config_dict['working mode']['relative error'] != 'None'):
                self.relative_error = float(config_dict['working mode']['relative error'])
                if (self.verbose >= 1):
                    self.logger.info('Stopping when relative error is below {0}'.format(self.relative_error))
            else:
                self.relative_error = 1e-4
        else:
            self.relative_error = 1e-4
        
        # Deal with the atmospheres
        tmp = config_dict['atmospheres']

        self.atmospheres = {}

        if (self.verbose >= 1):
            self.logger.info('Adding atmospheres')

        for key, value in tmp.items():
            
            if ('photosphere' in key):
                if (self.verbose >=1):
                    self.logger.info('  - New available photosphere : {0}'.format(value['name']))

                self.add_photosphere(value)
                                                            
            if ('chromosphere' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available chromosphere : {0}'.format(value['name']))

                self.add_chromosphere(value)
                                            
            if ('parametric' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available parametric : {0}'.format(value['name']))

                self.add_parametric(value)

            if ('straylight' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available straylight : {0}'.format(value['name']))

                self.add_straylight(value)

        self.setup()

    def setup(self):
        """
        Setup the model for synthesis/inversion. This setup includes adding the topologies, removing unused
        atmospheres, reading the number of cycles for the inversion and some sanity checks

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
                            
        # Adding topologies
        if (self.verbose >= 1):
            self.logger.info("Adding topologies") 
        for value in self.topologies:
            self.add_topology(value)

        # Remove unused atmospheres defined in the configuration file and not in the topology
        if (self.verbose >= 1):
            self.logger.info("Removing unused atmospheres")
        self.remove_unused_atmosphere()        

        # Calculate indices for atmospheres
        index_chromosphere = 1
        index_photosphere = 1
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                v.index = index_photosphere
                index_photosphere += 1
            if (v.type == 'chromosphere'):
                v.index = index_chromosphere
                index_chromosphere += 1

        # Check that number of pixels is the same for all atmospheric files if in synthesis mode
        if (self.working_mode == 'synthesis'):
            n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
            all_equal = all(x == n_pixels[0] for x in n_pixels)
            if (not all_equal):
                for k, v in self.atmospheres.items():
                    self.logger.info('{0} -> {1}'.format(k, v.n_pixel))
                raise Exception("Files with model atmospheres do not contain the same number of pixels")
            else:
                if (self.verbose >= 1):
                    self.logger.info('Number of pixels to read : {0}'.format(n_pixels[0]))
                self.n_pixels = n_pixels[0]

        if (self.working_mode == 'inversion'):
            n_pixels = [v.n_pixel for k, v in self.spectrum.items()]
            all_equal = all(x == n_pixels[0] for x in n_pixels)
            if (not all_equal):
                for k, v in self.spectrum.items():
                    self.logger.info('{0} -> {1}'.format(k, v.n_pixel))
                raise Exception("Files with spectral regions do not contain the same number of pixels")
            else:
                if (self.verbose >= 1):
                    self.logger.info('Number of pixels to invert : {0}'.format(n_pixels[0]))
                self.n_pixels = n_pixels[0]
        

        # Check that the number of pixels from all observations (in case of inversion) is the same
        # Check also if they are equal to those of the models
        # n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
        # all_equal = all(x == n_pixels[0] for x in n_pixels)

        # Check that the number of cycles is the same for all atmospheres (in case of inversion)
        if (self.working_mode == 'inversion'):
            cycles = []
            for k, v in self.atmospheres.items():
                for k2, v2 in v.cycles.items():                    
                    cycles.append(len(v2))

            all_equal = all(x == cycles[0] for x in cycles)
            if (not all_equal):
                raise Exception("Number of cycles in the nodes of active atmospheres is not always the same")
            else:
                if (self.n_cycles is None):
                    self.n_cycles = cycles[0]

        
        # if (self.working_mode == 'inversion'):
        #     cycles = []
        #     for tmp in ['I', 'Q', 'U', 'V']:
        #         if ( cycles.append
        #     for k, v in self.atmospheres.items():
        #         for k2, v2 in v.cycles.items():                    
        #             cycles.append(len(v2))

        #     all_equal = all(x == cycles[0] for x in cycles)
        #     if (not all_equal):
        #         raise Exception("Number of cycles in the nodes of active atmospheres is not always the same")
        #     else:
        #         if (self.n_cycles is None):
        #             self.n_cycles = cycles[0]

        filename = os.path.join(os.path.dirname(__file__),'data/LINEAS')
        ff = open(filename, 'r')
        self.LINES = ff.readlines()
        ff.close()

        self.init_sir()

        for k, v in self.spectrum.items():            
            v.allocate_info_cycles(n_cycles=self.n_cycles)

        for k, v in self.atmospheres.items():
            v.allocate_info_cycles(n_cycles=self.n_cycles)

        # Count total number of free parameters
        if (self.working_mode == 'inversion'):
            self.n_free_parameters = 0

            for k, v in self.atmospheres.items():
                for k2, v2 in v.cycles.items():                    
                    self.n_free_parameters += max(hazel.util.onlyint(v2[0:self.n_cycles+1]))

            if (self.verbose >= 1):
                self.logger.info('Total number of free parameters in all cycles : {0}'.format(self.n_free_parameters))
        
    def open_output(self):
        self.output_handler = Generic_output_file(self.output_file)        
        self.output_handler.open(self)

    def close_output(self):        
        self.output_handler.close()

    def write_output(self, randomization=0):
        if (self.working_mode == 'synthesis'):
            self.flatten_parameters_to_reference(cycle=0)        
        self.output_handler.write(self, pixel=0, randomization=randomization)

    def add_spectral(self, spectral):
        """
        Programmatically add a spectral region

        Parameters
        ----------
        spectral : dict
            Dictionary containing the following data
            'Name', 'Wavelength', 'Topology', 'Weights Stokes', 'Wavelength file', 'Wavelength weight file',
            'Observations file', 'Mask file'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        value = hazel.util.lower_dict_keys(spectral)

    
        if (self.verbose >= 1):            
            self.logger.info('Adding spectral region {0}'.format(value['name']))        

        if ('wavelength file' not in value):
            value['wavelength file'] = None
        elif (value['wavelength file'] == 'None'):
            value['wavelength file'] = None

        if ('wavelength weight file' not in value):
            value['wavelength weight file'] = None
        elif (value['wavelength weight file'] == 'None'):
            value['wavelength weight file'] = None

        if ('observations file' not in value):
            value['observations file'] = None
        elif (value['observations file'] == 'None'):
            value['observations file'] = None

        if ('stokes weights' not in value):
            value['stokes weights'] = None
        elif (value['stokes weights'] == 'None'):
            value['stokes weights'] = None

        if ('mask file' not in value):
            value['mask file'] = None
        elif (value['mask file'] == 'None'):
            value['mask file'] = None

        if ('los' not in value):
            value['los'] = None
        elif (value['los'] == 'None'):
            value['los'] = None

        for tmp in ['i', 'q', 'u', 'v']:
            if ('weights stokes {0}'.format(tmp) not in value):
                value['weights stokes {0}'.format(tmp)] = [None]*10
            elif (value['weights stokes {0}'.format(tmp)] == 'None'):
                value['weights stokes {0}'.format(tmp)] = [None]*10

        if ('boundary condition' not in value):
            value['boundary condition'] = None
        elif (value['boundary condition'] == 'None'):
            value['boundary condition'] = None

        if ('instrumental profile' not in value):
            value['instrumental profile'] = None
        elif (value['instrumental profile'] == 'None'):
            value['instrumental profile'] = None

        # Wavelength file is not present
        if (value['wavelength file'] is None):

            # If the wavelength is defined
            if ('wavelength' in value):
                axis = value['wavelength']
                wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
                if (self.verbose >= 1):
                    self.logger.info('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
            else:
                raise Exception('Wavelength range is not defined. Please, use "Wavelength" or "Wavelength file"')
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Reading wavelength axis from {0}'.format(value['wavelength file']))
            wvl = np.loadtxt(value['wavelength file'])
                
        if (value['wavelength weight file'] is None):
            if (self.verbose >= 1 and self.working_mode == 'inversion'):
                self.logger.info('  - Setting all wavelength weights to 1')
            weights = np.ones((4,len(wvl)))
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Reading wavelength weights from {0}'.format(value['wavelength weight file']))
            weights = np.loadtxt(value['wavelength weight file'], skiprows=1).T

        # Observations file not present
        if (value['observations file'] is None):
            if (self.working_mode == 'inversion'):
                raise Exception("Inversion mode without observations is not allowed.")            
            obs_file = None
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using observations from {0}'.format(value['observations file']))
            obs_file = value['observations file']

        if (value['mask file'] is None):            
            mask_file = None
            if (self.verbose >= 1):
                self.logger.info('  - No mask for pixels')
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using mask from {0}'.format(value['mask file']))
            mask_file = value['mask file']

        if (value['instrumental profile'] is None):
            if (self.verbose >= 1):
                self.logger.info('  - No instrumental profile')
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Instrumental profile : {0}'.format(value['instrumental profile']))

        # if (value['straylight file'] is None):
        #     if (self.verbose >= 1):
        #         self.logger.info('  - Not using straylight')
        #     stray_file = None
        # else:
        #     if (self.verbose >= 1):
        #         self.logger.info('  - Using straylight from {0}'.format(value['straylight file']))
        #     stray_file = value['straylight file']

        if (value['los'] is None):
            if (self.working_mode == 'synthesis'):
                raise Exception("You need to provide the LOS for spectral region {0}".format(value['name']))
            los = None
        else:
            los = np.array(value['los']).astype('float64')
            if (self.verbose >= 1):
                self.logger.info('  - Using LOS {0}'.format(value['los']))

        if (value['boundary condition'] is None):
            if (self.verbose >= 1):
                self.logger.info('  - Using default boundary conditions [1,0,0,0] in spectral region {0} or read from file. Check carefully!'.format(value['name']))
            boundary = np.array([1.0,0.0,0.0,0.0])            
            self.normalization = 'on-disk'
        else:
            boundary = np.array(value['boundary condition']).astype('float64')
            if (boundary[0] == 0.0):
                self.logger.info('  - Using off-limb normalization (peak intensity)')

            if (self.verbose >= 1):
                self.logger.info('  - Using boundary condition {0}'.format(value['boundary condition']))

        stokes_weights = []
        for st in ['i', 'q', 'u', 'v']:
            tmp = hazel.util.tofloat(value['weights stokes {0}'.format(st)])
            tmp = [i if i is not None else 1.0 for i in tmp]
            stokes_weights.append(tmp)
        
        stokes_weights = np.array(stokes_weights)

        self.spectrum[value['name']] = Spectrum(wvl=wvl, weights=weights, observed_file=obs_file, 
            name=value['name'], stokes_weights=stokes_weights, los=los, boundary=boundary, mask_file=mask_file, instrumental_profile=value['instrumental profile'])

        self.topologies.append(value['topology'])
        
    
    def add_photosphere(self, atmosphere):
        """
        Programmatically add a photosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = SIR_atmosphere(working_mode=self.working_mode, name=atm['name'])
        lines = [int(k) for k in list(atm['spectral lines'])]

        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        if ('reference frame' in atm):
            if ('line-of-sight' in atm['reference frame']):
                self.atmospheres[atm['name']].reference_frame = 'line-of-sight'
            if ('vertical' in atm['reference frame']):
                self.atmospheres[atm['name']].reference_frame = 'vertical'
        else:
            self.atmospheres[atm['name']].reference_frame = 'vertical'

        if (self.verbose >= 1):
            self.logger.info("    * Adding line : {0}".format(lines))
            self.logger.info("    * Magnetic field reference frame : {0}".format(self.atmospheres[atm['name']].reference_frame))
                    
        self.atmospheres[atm['name']].add_active_line(lines=lines, spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():            
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v                                            
        
        if ('reference atmospheric model' in atm):
            my_file = Path(atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()

        if ('nodes' in atm):            
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)        


    def add_chromosphere(self, atmosphere):
        """
        Programmatically add a chromosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)
        
        self.atmospheres[atm['name']] = Hazel_atmosphere(working_mode=self.working_mode, name=atm['name'])

        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        self.atmospheres[atm['name']].add_active_line(line=atm['line'], spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        if ('reference frame' in atm):
            if (atm['reference frame'] == 'line-of-sight'):
                self.atmospheres[atm['name']].reference_frame = 'line-of-sight'
            if (atm['reference frame'] == 'vertical'):
                self.atmospheres[atm['name']].reference_frame = 'vertical'
        else:
            self.atmospheres[atm['name']].reference_frame = 'vertical'

        if (self.verbose >= 1):
            self.logger.info("    * Adding line : {0}".format(atm['line']))
            self.logger.info("    * Magnetic field reference frame : {0}".format(self.atmospheres[atm['name']].reference_frame))

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v
        
        if ('coordinates for magnetic field vector' in atm):
            if (atm['coordinates for magnetic field vector'] == 'cartesian'):
                self.atmospheres[atm['name']].coordinates_B = 'cartesian'
            if (atm['coordinates for magnetic field vector'] == 'spherical'):
                self.atmospheres[atm['name']].coordinates_B = 'spherical'
        else:
            self.atmospheres[atm['name']].coordinates_B = 'cartesian'

        self.atmospheres[atm['name']].select_coordinate_system()

        if (self.verbose >= 1):            
            self.logger.info("    * Magnetic field coordinates system : {0}".format(self.atmospheres[atm['name']].coordinates_B))            


        if ('reference atmospheric model' in atm):
            my_file = Path(atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        self.atmospheres[atm['name']].height = float(atm['height'])

        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)
                            
            
    def add_parametric(self, atmosphere):
        """
        Programmatically add a parametric atmosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Wavelength', 'Reference atmospheric model', 'Type',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Parametric_atmosphere(working_mode=self.working_mode)

        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        if ('reference atmospheric model' in atm):
            my_file = Path(atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v

    
    def add_straylight(self, atmosphere):
        """
        Programmatically add a straylight atmosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Straylight_atmosphere(working_mode=self.working_mode)
                
        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)        

        my_file = Path(atm['reference atmospheric model'])
        if (not my_file.exists()):
            raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

        if ('reference atmospheric model' in atm):
            self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v
        

    def remove_unused_atmosphere(self):
        """
        Remove unused atmospheres
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        
        to_remove = []        
        for k, v in self.atmospheres.items():
            if (not v.active):
                to_remove.append(k)
                if (self.verbose >= 1):
                    self.logger.info('  - Atmosphere {0} deleted.'.format(k))
                
        for k in to_remove:
            self.atmospheres.pop(k)
                    
    def init_sir_external(self):
        """
        Initialize SIR for this synthesis
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                f = open('lte.grid', 'w')
                f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
                f.write("           b) The first six characters of the last line                \n")
                f.write("          in the header (if any) must contain the symbol ---       \n")
                f.write("\n")                                                                       
                f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
                f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
                f.write("-----------------------------------------------------------------------\n")

                ind_low = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[0])).argmin()
                ind_top = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[1])).argmin()

                low = v.spectrum.wavelength_axis[ind_low]
                top = v.spectrum.wavelength_axis[ind_top]         # TODO
                delta = (v.spectrum.wavelength_axis[1] - v.spectrum.wavelength_axis[0])

                filename = os.path.join(os.path.dirname(__file__),'data/LINEAS')
                ff = open(filename, 'r')
                flines = ff.readlines()
                ff.close()

                for i in range(len(v.lines)):
                    for l in flines:
                        tmp = l.split()
                        index = int(tmp[0].split('=')[0])
                        if (index == v.lines[0]):
                            wvl = float(tmp[2])                    
                                            
                f.write("{0}            :  {1}, {2}, {3}\n".format(str(v.lines)[1:-1], 1e3*(low-wvl), 1e3*delta, 1e3*(top-wvl)))
                f.close()
                
                v.n_lambda = sir_code.init_externalfile(v.index, filename)

    def init_sir(self):
        """
        Initialize SIR for this synthesis. This version does not make use of any external file, which might be
        not safe when running in MPI mode.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        lines = []
        n_lines = 0

        elements = {'H':1,'HE':2,'LI':3,'BE':4,'B':5,'C':6,'N':7,'O':8,'F':9,'NE':10,
            'NA':11,'MG':12,'AL':13,'SI':14,'P':15,'S':16,'CL':17,'AR':18,'K':19,'CA':20,'SC':21,'TI':22,'V':23,'CR':24,
            'MN':25,'FE':26,'CO':27,'NI':28,'CU':29,'ZN':30,'GA':31,'GE':32,'AS':33,'SE':34,'BR':35,'KR':36,
            'RB':37,'SR':38,'Y':39,'ZR':40,'NB':41,'MO':42,'TC':43,'RU':44,'RH':45,'PD':46,'AG':47,'CD':48,'IN':49,
            'SN':50,'SB':51,'TE':52,'I':53,'XE':54,'CS':55,'BA':56,'LA':57,'CE':58,'PR':59,'ND':60,'PM':61,
            'SM':62,'EU':63,'GD':64,'TB':65,'DY':66,'HO':67,'ER':68,'TM':69,'YB':70,'LU':71,'HF':72,'TA':73,'W':74,
            'RE':75,'OS':76,'IR':77,'PT':78,'AU':79,'HG':80,'TL':81,'PB':82,'BI':83,'PO':84,'AT':85,'RN':86,
            'FR':87,'RA':88,'AC':89,'TH':90,'PA':91,'U':92}
        states = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}

        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):

                n_lines += 1
                
                ind_low = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[0])).argmin()
                ind_top = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[1])).argmin()

                low = v.spectrum.wavelength_axis[ind_low]
                top = v.spectrum.wavelength_axis[ind_top]         # TODO
                delta = (v.spectrum.wavelength_axis[1] - v.spectrum.wavelength_axis[0])
                
                nblend = len(v.lines)

                lines = np.zeros(len(v.lines), dtype=np.intc)
                atom = np.zeros(len(v.lines), dtype=np.intc)
                istage = np.zeros(len(v.lines), dtype=np.intc)
                wvl = np.zeros(len(v.lines))
                zeff = np.zeros(len(v.lines))
                energy = np.zeros(len(v.lines))
                loggf = np.zeros(len(v.lines))
                mult1 = np.zeros(len(v.lines), dtype=np.intc)
                mult2 = np.zeros(len(v.lines), dtype=np.intc)
                design1 = np.zeros(len(v.lines), dtype=np.intc)
                design2 = np.zeros(len(v.lines), dtype=np.intc)
                tam1 = np.zeros(len(v.lines))
                tam2 = np.zeros(len(v.lines))
                alfa = np.zeros(len(v.lines))
                sigma = np.zeros(len(v.lines))
                
                for i in range(len(v.lines)):            
                    lines[i] = v.lines[i]                    
                    for l in self.LINES:
                        tmp = l.split()
                        index = int(tmp[0].split('=')[0])
                        if (index == v.lines[i]):
                                                        
                            atom[i] = elements[tmp[0].split('=')[1]]
                            istage[i] = tmp[1]
                            wvl[i] = float(tmp[2])
                            zeff[i] = float(tmp[3])
                            energy[i] = float(tmp[4])
                            loggf[i] = float(tmp[5])
                            mult1[i] = int(tmp[6][:-1])
                            mult2[i] = int(tmp[8][:-1])
                            design1[i] = states[tmp[6][-1]]
                            design2[i] = states[tmp[8][-1]]
                            tam1[i] = float(tmp[7].split('-')[0])
                            tam2[i] = float(tmp[9].split('-')[0])
                            if (len(tmp) == 12):
                                alfa[i] = float(tmp[-2])
                                sigma[i] = float(tmp[-1])
                            else:
                                alfa[i] = 0.0
                                sigma[i] = 0.0
                
                lambda0 = 1e3*(low-wvl[0])
                lambda1 = 1e3*(top-wvl[0])
                n_steps = ind_top - ind_low + 1

                v.n_lambda = n_steps
                
                sir_code.init(v.index, nblend, lines, atom, istage, wvl, zeff, energy, loggf,
                    mult1, mult2, design1, design2, tam1, tam2, alfa, sigma, lambda0, lambda1, n_steps)

    def exit_hazel(self):
        for k, v in self.atmospheres.items():            
            if (v.type == 'chromosphere'):
                hazel_code.exit(v.index)

    def add_topology(self, atmosphere_order):
        """
        Add a new topology

        Parameters
        ----------
        topology : str
            Topology
        
        Returns
        -------
        None

        """

        # Transform the order to a list of lists
        if (self.verbose >= 1):
            self.logger.info('  - {0}'.format(atmosphere_order))

        vertical_order = atmosphere_order.split('->')        
        order = []
        for k in vertical_order:
            name = k.strip().replace('(','').replace(')','').split('+')
            name = [k.strip() for k in name]
            
            tmp = []
            for n in name:
                tmp.append(n)
                self.atmospheres[n].active = True

            order.append(tmp)
        
        order_flat = [item for sublist in order for item in sublist]

        # Check that straylight components, if any, are not at the last position
        for atm in order_flat[:-1]:            
            if (self.atmospheres[atm].type == 'straylight'):
                raise Exception("Straylight components can only be at the last position of a topology.")
        
        self.order_atmospheres.append(order)
                
    def normalize_ff(self):
        """
        Normalize all filling factors so that they add to one to avoid later problems.
        We use a softmax function to make sure they all add to one and can be unconstrained

        ff_i = exp(x_i) / sum(exp(x_i))

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """

        for atmospheres in self.order_atmospheres:
            for order in atmospheres:

                total_ff = 0.0
                for atm in order:            
                    if (self.atmospheres[atm].type is not 'straylight'):
                        if (self.working_mode == 'inversion'):                            
                            self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1]
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                        else:
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)
                        total_ff += ff

                for atm in order:                    
                    if (self.atmospheres[atm].type is not 'straylight'):
                        if (self.working_mode == 'inversion'):
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                            self.atmospheres[atm].parameters['ff'] = ff / total_ff
                            self.atmospheres[atm].parameters['ff'] = physical_to_transformed(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                        else:
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)
                            self.atmospheres[atm].parameters['ff'] = ff / total_ff
                            self.atmospheres[atm].parameters['ff'] = physical_to_transformed(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)


    def synthesize_spectral_region(self, spectral_region, perturbation=False):
        """
        Synthesize all atmospheres for a single spectral region and normalize to the continuum of the quiet Sun at disk center

        Parameters
        ----------
        spectral_region : str
            Spectral region to synthesize
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """
        stokes = None
        stokes_out = None

        # Loop over all atmospheres        
        for i, atmospheres in enumerate(self.order_atmospheres):
            
            for n, order in enumerate(atmospheres):
                                                                
                for k, atm in enumerate(order):                    
                    
                    if (self.atmospheres[atm].spectrum.name == spectral_region):                        
                        
                        # Update the boundary condition only for the first atmosphere if several are sharing ff      
                        if (n > 0 and k == 0):
                            ind_low, ind_top = self.atmospheres[atm].wvl_range
                            
                            if (perturbation):
                                stokes_out = self.atmospheres[atm].spectrum.stokes_perturbed[:, ind_low:ind_top] * hazel.util.i0_allen(self.atmospheres[atm].spectrum.wavelength_axis[ind_low:ind_top], 1.0)[None,:]
                            else:
                                stokes_out = self.atmospheres[atm].spectrum.stokes[:, ind_low:ind_top] * hazel.util.i0_allen(self.atmospheres[atm].spectrum.wavelength_axis[ind_low:ind_top], 1.0)[None,:]

                        if (self.atmospheres[atm].type == 'straylight'):
                            stokes, error = self.atmospheres[atm].synthesize()
                            if (error == 1):
                                raise 
                            stokes += (1.0 - self.atmospheres[atm].parameters['ff']) * stokes_out                            
                        else:
                            if (k == 0):                                
                                stokes, error = self.atmospheres[atm].synthesize(stokes_out)
                            else:             
                                tmp, error = self.atmospheres[atm].synthesize(stokes_out)       
                                stokes += tmp
                                    
                        ind_low, ind_top = self.atmospheres[atm].wvl_range
                        
                        if (perturbation):
                            self.atmospheres[atm].spectrum.stokes_perturbed[:, ind_low:ind_top] = stokes / hazel.util.i0_allen(self.atmospheres[atm].spectrum.wavelength_axis[ind_low:ind_top], 1.0)[None,:]
                        else:                            
                            self.atmospheres[atm].spectrum.stokes[:, ind_low:ind_top] = stokes / hazel.util.i0_allen(self.atmospheres[atm].spectrum.wavelength_axis[ind_low:ind_top], 1.0)[None,:]
    
    def synthesize(self, perturbation=False):
        """
        Synthesize all atmospheres

        Parameters
        ----------
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """
        if (self.working_mode == 'inversion'):
            self.normalize_ff()
        for k, v in self.spectrum.items():            
            self.synthesize_spectral_region(k, perturbation=perturbation)
            if (v.normalization == 'off-limb'):
                if (perturbation):
                    v.stokes_perturbed /= np.max(v.stokes_perturbed[0,:])
                else:
                    v.stokes /= np.max(v.stokes[0,:])

            if (v.psf_spectral is not None):
                for i in range(4):                    
                    v.stokes[i,:] = scipy.signal.convolve(v.stokes[i,:], v.psf_spectral, mode='same', method='auto')
                
            
    def find_active_parameters(self, cycle):
        """
        Find all active parameters in all active atmospheres in the current cycle

        Parameters
        ----------
        cycle : int
            Cycle to consider
        
        Returns
        -------
        None

        """
        pars = []
        coupled = []
        self.nodes = []
        left = 0
        right = 0
        for atmospheres in self.order_atmospheres:
            for n, order in enumerate(atmospheres):
                for k, atm in enumerate(order):                    
                    for l, par in self.atmospheres[atm].cycles.items():
                        if (par is not None):
                            if (hazel.util.isint(par[cycle])):
                                if (par[cycle] > 0):
                                    
                                    # [Atmosphere name, n_nodes, nodes, value, range]
                                    self.atmospheres[atm].nodes[l] = np.zeros(par[cycle])

                                    self.atmospheres[atm].n_nodes[l] = par[cycle]

                                    right += par[cycle]
                                    
                                    n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                    tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 
                                        'ranges': self.atmospheres[atm].ranges[l], 'delta': self.atmospheres[atm].epsilon[l],
                                        'left': left, 'right': right, 'regularization': self.atmospheres[atm].regularization[l],
                                        'coupled': False}

                                    self.nodes.append(self.atmospheres[atm].nodes[l])
                                        
                                    left = copy.copy(right)
                                                                                                            
                                    pars.append(tmp)

                                else:
                                    self.atmospheres[atm].nodes[l] = 0.0
                                    self.atmospheres[atm].n_nodes[l] = 0
                            else:
                                
                                n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 'coupled': True}

                                coupled.append(tmp)

        self.active_meta = pars
        self.coupled_meta = coupled

        if (not self.nodes):
            raise Exception("No parameters to invert in cycle {0}. Please add them or reduce the number of cycles. ".format(cycle))            

        self.nodes = np.concatenate(self.nodes).ravel()
        
        
    def synthesize_and_compute_rf(self, compute_rf=False, include_jacobian=False):
        """
        Compute response functions for all free parameters according to all active_parameters

        Parameters
        ----------
        compute_rf : bool (optional, default False)
            If True, then compute the response functions. If not, just compute the synthesis.
        
        Returns
        -------
        None

        """

        self.synthesize()

        if (not compute_rf):            
            return

        n_active_pars = len(self.active_meta)

        loop = 0

        self.hessian_regularization = np.zeros(self.n_free_parameters_cycle)
        self.grad_regularization = np.zeros(self.n_free_parameters_cycle)

        for par in self.active_meta:            
            nodes = self.nodes[par['left']:par['right']]

            lower = par['ranges'][0]
            upper = par['ranges'][1]

            if (self.verbose >= 4):
                self.logger.info(" * RF to {0} - {1} - nodes={2}".format(par['parameter'], par['atm'], par['n_nodes']))
                        
            for i in range(par['n_nodes']):
                perturbation = np.zeros(par['n_nodes'])
                if (nodes[i] == 0):
                    perturbation[i] = self.epsilon * par['delta']
                else:
                    perturbation[i] = self.epsilon * nodes[i]                

                # Perturb this parameter
                self.atmospheres[par['atm']].nodes[par['parameter']] = nodes + perturbation

                # Also perturb those parameters that are coupled
                for par2 in self.coupled_meta:
                    if (par2['coupled'] is True):                            
                        if (par['atm'] == par2['n_nodes'] and par['parameter'] == par2['parameter']):
                            if (self.verbose >= 4):
                                self.logger.info("   * Coupling RF to {0} - {1}".format(par2['parameter'], par2['atm']))
                            self.atmospheres[par2['atm']].nodes[par2['parameter']] = nodes + perturbation
                                
                # Synthesize
                self.synthesize(perturbation=True)
                                                
                # And come back to the original value of the nodes
                self.atmospheres[par['atm']].nodes[par['parameter']] = nodes

                for par2 in self.coupled_meta:
                    if (par2['coupled'] is True):
                        if (par['atm'] == par2['n_nodes'] and par['parameter'] == par2['parameter']):
                            self.atmospheres[par2['atm']].nodes[par2['parameter']] = nodes
                
                if (include_jacobian):
                    # jacobian =
                    # self.atmospheres[par['atm']].jacobian[par['parameter']]                    
                    jacobian = jacobian_transformation(nodes[i], lower, upper)                    
                else:
                    jacobian = 1.0

                rf = {}
                for k, v in self.spectrum.items():
                    rf[k] = jacobian * np.expand_dims((v.stokes - v.stokes_perturbed) / perturbation[i], 0)

                # rf = np.expand_dims((self.spectrum['spec1'].stokes - self.spectrum['spec1'].stokes_perturbed) / perturbation[i], 0)
                                                
                if (loop == 0):
                    self.response = rf
                else:
                    # self.response = np.vstack([self.response, rf])
                    for k, v in self.spectrum.items():
                        self.response[k] = np.vstack([self.response[k], rf[k]])

                if (par['regularization'] is not None):
                    if (par['regularization'][0] == 'l2-value'):
                        alpha = float(par['regularization'][1])
                        lower = par['ranges'][0]
                        upper = par['ranges'][1]
                        value = physical_to_transformed(float(par['regularization'][2]), lower, upper)
                        self.grad_regularization[par['left']:par['right']] = 2.0 * alpha * (self.atmospheres[par['atm']].nodes[par['parameter']] - value)
                        self.hessian_regularization[par['left']:par['right']] = 2.0 * alpha

                loop += 1

        
    def flatten_parameters_to_reference(self, cycle):
        """
        Flatten all current parameters to the reference atmosphere

        Parameters
        ----------
        cycle : int
            Current cycle
        
        Returns
        -------
        None

        """                
        if (self.working_mode == 'inversion'):
            for k, v in self.atmospheres.items():
                v.set_reference(cycle=cycle)
            
        for k, v in self.spectrum.items():            
            v.stokes_cycle[cycle] = copy.deepcopy(v.stokes)

            if (self.working_mode == 'inversion'):                
                v.chi2_cycle[cycle] = copy.deepcopy(v.chi2)                        
                v.bic_cycle[cycle] = copy.deepcopy(self.n_free_parameters * np.log(v.dof) + v.dof * np.log(v.rss))
                v.aic_cycle[cycle] = copy.deepcopy(2.0 * self.n_free_parameters + v.dof * np.log(v.rss))
        
    def set_new_model(self, nodes):
        """
        Set the nodes of the current model to the values passed on the arguments

        Parameters
        ----------
        nodes : float
            Array with the new set of nodes
        
        Returns
        -------
        None

        """

        n_active_pars = len(self.active_meta)  

        # Modify all active parameters
        for par in self.active_meta:
            left = par['left']
            right = par['right']

            self.atmospheres[par['atm']].nodes[par['parameter']] = nodes[left:right]

        # Modify all coupled parameters accordingly
        for par in self.coupled_meta:
            for par2 in self.active_meta:
                if (par2['atm'] == par['n_nodes'] and par2['parameter'] == par['parameter']):
                    
                    left = par2['left']
                    right = par2['right']
                    
                    self.atmospheres[par['atm']].nodes[par['parameter']] = nodes[left:right]
                    self.atmospheres[par['atm']].parameters[par['parameter']] = copy.copy(self.atmospheres[par2['atm']].parameters[par2['parameter']])
                    

    def modified_svd_inverse(self, H, tol=1e-8):
        """
        Compute the inverse of the Hessian matrix using a modified SVD, by thresholding each subpsace separately

        Parameters
        ----------
        H : float
            Hessian matrix

        tol : float
            Tolerance for the singular value of each subspace
        
        Returns
        -------
        None

        """

        try:
            U, w, VT = np.linalg.svd(H, full_matrices=False)
        except np.linalg.LinAlgError:
            U, w, VT = scipy.linalg.svd(H, full_matrices=False, lapack_driver='gesvd')   # This calculation should be more robust but slower

        w_new = np.zeros_like(w)
        
        for par in self.active_meta:            
            left = par['left']
            right = par['right']

            Ui = np.zeros_like(U)
            Ui[:,left:right] = U[:,left:right]

            Gamma_i = np.diagonal(np.diag(w) @ Ui.T @ U).copy()
            
            wmax = np.max(np.abs(Gamma_i))            
            Gamma_i[np.abs(Gamma_i) < tol*wmax] = 0.0

            w_new += Gamma_i

        w_new_inv = np.zeros_like(w)
        ind = np.where(w_new != 0)[0]
        w_new_inv[ind] = 1.0 / w_new[ind]
        
        return U, w_new_inv, VT
        

    def compute_chi2(self, only_chi2=False, weights=None):
        """
        Compute chi2 for all spectral regions

        Parameters
        ----------
        obs : float
            Vector of observations
        only_chi2 : bool
            Control whether the gradient and Hessian is returned
        
        Returns
        -------
        None

        """
        chi2 = 0.0
        rss = 0.0
        n = len(self.nodes)
        dchi2 = np.zeros(n)
        ddchi2 = np.zeros((n,n))

        for k, v in self.spectrum.items():
            residual = (v.stokes - v.obs)
            
            # Do not use weights. This is used for the computation of errors            
            # if (weights is None):
            weights = (v.stokes_weights[:,self.cycle][:,None] * v.wavelength_weights) * v.factor_chi2            

            chi2 += np.sum(weights * residual**2)

            rss += np.sum(residual**2)
            
            if (not only_chi2):
                response = self.response[k]
                dchi2 += -2.0 * np.sum(weights[None,:,:] * response * residual[None,:,:] , axis=(1,2)) #/ v.dof
                ddchi2 += 2.0 * np.sum(weights[None,None,:,:] * response[None,:,:,:] * response[:,None,:,:] , axis=(2,3)) #/ v.dof

            v.chi2 = chi2
            v.rss = rss
            
        if (not only_chi2):                
            return chi2, dchi2, ddchi2
        else:                
            return chi2

        # if (not only_chi2):            
        #     return chi2, dchi2, ddchi2
        # else:
        #     return chi2

    def compute_uncertainty(self):
        """
        Compute the uncertainty in the parameters at the minimum with the current Hessian

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
    
        #----------------------------
        # Recalculate chi2 without weights
        # chi2 = 0.0
        # for k, v in self.spectrum.items():
        #     residual = (v.stokes - v.obs)

        # weights = v.dof / residual**2

        # # Calculate Hessian
        # self.synthesize_and_compute_rf(compute_rf=True)
        # chi2, dchi2, ddchi2 = self.compute_chi2(weights=weights)
        # hessian = 0.5 * ddchi2


        #----------------------------
        # Recalculate chi2 without weights        

        # Calculate Hessian
        self.synthesize_and_compute_rf(compute_rf=True, include_jacobian=True)
        chi2, dchi2, ddchi2 = self.compute_chi2()
        hessian = 0.5 * ddchi2
        
        U, w_inv, VT = self.modified_svd_inverse(hessian, tol=1e-12)
        cov = VT.T.dot(np.diag(w_inv)).dot(U.T)

        # breakpoint()

        for par in self.active_meta:
            left = par['left']
            right = par['right']

            dof = self.atmospheres[par['atm']].spectrum.dof
            rf = scipy.stats.chi2(dof)
            delta = np.sqrt(rf.isf(1.0 - scipy.special.erf(1.0/np.sqrt(2.0))))
            
            rf = scipy.stats.chi2(right-left)
            delta = np.sqrt(rf.isf(1.0 - scipy.special.erf(1.0/np.sqrt(2.0))))

            cov_diagonal = np.abs(np.diagonal(cov[left:right,left:right]))

            # This gives 1sigma error in the transformed domain
            error = np.sqrt(cov_diagonal) * delta

            # Multiply by the Jacobian of the transformation to compute the error in the physical quantities
            error *= jacobian_transformation(self.nodes[left:right], par['ranges'][0], par['ranges'][1])
                    
            self.atmospheres[par['atm']].error[par['parameter']] = error

    def _fun_backtracking(self, log_lambda, dchi2, ddchi2):
        H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
        H += np.diag(10.0**(log_lambda) * np.diag(H))
        gradF = 0.5 * (dchi2 + self.grad_regularization)

        U, w_inv, VT = self.modified_svd_inverse(H, tol=1e-8)

        # xnew = xold - H^-1 * grad F
        delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
        
        # Clip the new solution so that the step is resaonable
        new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
    
        self.set_new_model(new_solution)
                    
        self.synthesize_and_compute_rf()
        
        chi2 = self.compute_chi2(only_chi2=True)

        if (self.verbose >= 4):
            self.logger.info('  - Backtracking - lambda: {0:7.5f} - chi2: {1:7.5f}'.format(10.0**log_lambda, chi2))
                    
        return chi2
            
    
    def backtracking_brent(self, dchi2, ddchi2, maxiter=10, bounds=[-3.0,3.0], tol=1e-2):

        tmp = scipy.optimize.minimize_scalar(self._fun_backtracking, bounds=bounds, args=(dchi2, ddchi2), 
            method='bounded', tol=tol, options={'maxiter': maxiter})

        return 10.0**tmp['x']
        
    def backtracking_parabolic(self, dchi2, ddchi2, direction='down', maxiter=5, lambda_init=1e-3, current_chi2=1e10):
        """
        Do the backtracking to get an optimal value of lambda in the LM algorithm

        Parameters
        ----------
        dchi2 : float
            Gradient of the chi2
        ddchi2 : float
            Second order derivatives with which the Hessian is computed
        direction : str, optional
            Direction on which do the backtracking ('down'/'up' for decreasing/increasing lambda)
        maxiter : int
            Maximum number of iterations
        lambda_init : float
            Initial value of lambda
        current_chi2 : float
            Current best chi2 to compare with those of the backtracking
        
        Returns
        -------
        lambda_opt : float
            Optimal value of lambda found. Bracketed value if bracketing has been possible or just the best value otherwise
        bracketed : bool
            True if the best value has been bracketed
        best_chi2 : float
            Best value of chi2 found
        """
        
        lambdaLM = lambda_init

        chi2_arr = []
        lambdas = []
        sols = []
        keepon = True
        bracketed = False
        loop = 0
        best_chi2 = current_chi2

        while keepon:

            H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
            H += np.diag(lambdaLM * np.diag(H))
            gradF = 0.5 * (dchi2 + self.grad_regularization)

            U, w_inv, VT = self.modified_svd_inverse(H, tol=1e-8)

            # xnew = xold - H^-1 * grad F
            delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
            
            # Clip the new solution so that the step is resaonable
            new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
            sols.append(new_solution)
        
            self.set_new_model(new_solution)
                        
            self.synthesize_and_compute_rf()
            
            chi2_arr.append(self.compute_chi2(only_chi2=True))
            
            lambdas.append(lambdaLM)

            if (self.verbose >= 4):
                if (direction == 'down'):
                    self.logger.info('  - Backtracking: {0:2d} - lambda: {1:7.5f} - chi2: {2:7.5f}'.format(loop, lambdaLM, chi2_arr[-1]))
                else:
                    self.logger.info('  * Backtracking: {0:2d} - lambda: {1:7.5f} - chi2: {2:7.5f}'.format(loop, lambdaLM, chi2_arr[-1]))
            
            # If we improve the chi2
            if (chi2_arr[-1] < best_chi2):
                best_chi2 = chi2_arr[-1]

            ind_min = np.argmin(chi2_arr)
            
            if (loop > 1):
                
                # Have we bracketed the minimum
                if (ind_min != 0 and ind_min != len(chi2_arr)-1):
                    keepon = False
                    bracketed = True

            # If lambda < 1e-3, then stop
            if (lambdaLM < 1e-3 or loop > maxiter):
                keepon = False
                min_found = False
            
            if (direction == 'down'):
                lambdaLM /= np.sqrt(10.0)
            else:
                lambdaLM *= np.sqrt(10.0)
            
            loop += 1

        # Parabolic interpolation of the optimal value of lambda
        if (bracketed):
            coeff = np.polyfit(np.log(lambdas[ind_min-1:ind_min+2]), chi2_arr[ind_min-1:ind_min+2], 2)
            lambda_opt = np.exp(-coeff[1] / (2.0*coeff[0]))
        else:
            lambda_opt = lambdas[ind_min]

        return lambda_opt, bracketed, best_chi2, np.min(chi2_arr)

    def randomize(self):
        """
        Randomize all free parameters to lie uniformly in the interval [-2,2] in the transformed
        domain
        """        
        self.nodes = np.random.uniform(low=-2.0, high=2.0, size=self.nodes.shape)

    def invert(self, randomize=False, randomization_ind=None):
        """
        Invert all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        
        first = True


        # Reset reference model to the one loaded from the file
        for k, v in self.atmospheres.items():
            v.reset_reference()

        # Compute normalization factor for the chi^2    
        for k, v in self.spectrum.items():
            v.factor_chi2 = 1.0 / (v.noise**2 * v.dof)        
        
        lambdaLM = 10.0
        lambda_opt = 10.0
        bestchi2 = 1e10
    
        for self.cycle in range(self.n_cycles):            

            if (self.verbose >= 2):
                self.logger.info('-------------')
                if (randomization_ind):
                    self.logger.info('  Cycle {0} - Randomization {1} '.format(self.cycle, randomization_ind))
                else:
                    self.logger.info('  Cycle {0}  '.format(self.cycle))
                
                for k, v in self.spectrum.items():
                    self.logger.info('  Weights for region {0} : SI={1} - SQ={2} - SU={3} - SV={4}'.format(k, v.stokes_weights[0,self.cycle], v.stokes_weights[1,self.cycle],
                        v.stokes_weights[2,self.cycle], v.stokes_weights[3,self.cycle]))
                                

                self.logger.info('-------------')
                

            # Find all active parameters for this cycle and print them in the output
            self.find_active_parameters(self.cycle)

            tmp = [pars['atm'] for pars in self.active_meta]
            tmp = list(set(tmp))

            self.n_free_parameters_cycle = 0

            
            for k, v in self.atmospheres.items():
                if (k in tmp):
                    if (self.verbose >= 3):
                        self.logger.info('Free parameters for {0}'.format(k))
                    for pars in self.active_meta:
                        if (pars['atm'] == k):
                            if (self.verbose >= 3):                                
                                if (pars['coupled'] is False):
                                    if (pars['n_nodes'] == 1):
                                        if (pars['regularization'] is not None):
                                            self.logger.info('  - {0} with {1} node - Regularization -> type:{2}, weight:{3}, value:{4}'.format(pars['parameter'], 
                                                pars['n_nodes'], pars['regularization'][0], pars['regularization'][1], pars['regularization'][2]))
                                        else:
                                            self.logger.info('  - {0} with {1} node - Not regularized'.format(pars['parameter'], pars['n_nodes']))
                                    else:
                                        if (pars['regularization'] is not None):
                                            self.logger.info('  - {0} with {1} nodes - Regularization -> type:{2}, weight:{3}, value:{4}'.format(pars['parameter'], 
                                                pars['n_nodes'], pars['regularization'][0], pars['regularization'][1], pars['regularization'][2]))
                                        else:
                                            self.logger.info('  - {0} with {1} nodes - Not regularized'.format(pars['parameter'], pars['n_nodes']))
                                else:
                                    self.logger.info('  - {0} coupled to {1} variable'.format(pars['parameter'], pars['n_nodes']))
                            if (pars['coupled'] is False):
                                self.n_free_parameters_cycle += pars['n_nodes']
            
            # Randomize parameters if necessary
            if (randomize):
                self.randomize()

            keepon = True
            iteration = 0

            # Mail Levenberg-Marquardt algorithm
            self.synthesize_and_compute_rf(compute_rf=True)
            chi2, dchi2, ddchi2 = self.compute_chi2()

            while keepon:                                
                
                # Simple parabolic backtracking
                if (self.backtracking == 'parabolic'):
                    lambda_opt, bracketed, best_chi2, backtracking_bestchi2_down = self.backtracking(dchi2, ddchi2, direction='down', maxiter=5, lambda_init=lambdaLM, current_chi2=chi2)

                    backtracking_bestchi2 = copy.copy(backtracking_bestchi2_down)

                    # If solution is not bracketed, then try on the other sense and use the best of the two
                    if (not bracketed):
                        lambda_opt_up, bracketed, best_chi2_up, backtracking_bestchi2_up = self.backtracking(dchi2, ddchi2, direction='up', maxiter=2, lambda_init=lambdaLM)                    

                        if (best_chi2_up < best_chi2):
                            lambda_opt = lambda_opt_up

                        backtracking_bestchi2 = np.min([backtracking_bestchi2, backtracking_bestchi2_up])

                # Bounded Brent backtracking
                if (self.backtracking == 'brent'):                    
                    lambda_opt = self.backtracking_brent(dchi2, ddchi2, maxiter=10, bounds=[-4.0,4.0], tol=1e-2)
                                                
                # if (self.verbose >= 3):
                    # self.logger.info('  * Optimal lambda: {0}'.format(lambda_opt))
                                            
                # If after backtracking the chi2 is larger than the current one, then increase lambda and go to the iteration
                # print(chi2, backtracking_bestchi2)
                # if (chi2 < backtracking_bestchi2 and iteration > 1):
                #     lambdaLM *= 100.0
                #     # print('breaking')
                #     continue


                # Give the final step
                H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
                H += np.diag(lambda_opt * np.diag(H))
                gradF = 0.5 * (dchi2 + self.grad_regularization)

                U, w_inv, VT = self.modified_svd_inverse(H, tol=1e-8)

                # xnew = xold - H^-1 * grad F
                delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)

                # New solution
                # Clip the new solution so that the step is resaonable
                new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
                self.set_new_model(new_solution)

                # Clip the new solution so that the step is resaonable
                self.nodes += np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)

                self.synthesize_and_compute_rf(compute_rf=True)
                chi2, dchi2, ddchi2 = self.compute_chi2()

                rel = 2.0 * (chi2 - bestchi2) / (chi2 + bestchi2)

                if (self.verbose > 2):
                    for k, v in self.atmospheres.items():
                        self.logger.info('')
                        self.logger.info('-----------')
                        self.logger.info('{0}'.format(k))
                        self.logger.info('-----------')
                        if (v.type == 'chromosphere'):
                            v.print_parameters(first=first)
                        if (v.type == 'photosphere'):
                            v.print_parameters(first=first)                        
                    first = False

                if (self.verbose >= 2):
                    self.logger.info('==============================================================================')
                    self.logger.info('It: {0} - chi2: {1:10.6f} - lambda_opt: {2:10.6f} - rel: {3:10.6f}'.format(iteration, chi2, lambda_opt, np.abs(rel)))
                    self.logger.info('==============================================================================')

                # Increase the optimal by 100 to find again the optimal value
                lambdaLM = 100.0 * lambda_opt

                bestchi2 = copy.copy(chi2)

                if (np.abs(rel) < self.relative_error or iteration > self.max_iterations):
                    keepon = False

                iteration += 1

                
                        
            self.set_new_model(self.nodes)

            # Calculate final chi2
            # self.synthesize_and_compute_rf()
            # chi2 = self.compute_chi2(only_chi2=True)
            

            self.compute_uncertainty()
            # if (self.verbose >= 2):
            #     self.atmospheres['ch1'].print_parameters(first=first, error=True)

            self.flatten_parameters_to_reference(self.cycle)

    def _func_grad(self, x):
        """
        Auxiliary functions to use with optimization methods that use gradients
        """
        self.nodes = x
        self.set_new_model(self.nodes)
        self.synthesize_and_compute_rf(compute_rf=True)
        self.chi2, dchi2, _ = self.compute_chi2()
        return self.chi2, dchi2

    def _func_nograd(self, x):
        """
        Auxiliary functions to use with optimization methods that do not use gradients
        """        
        self.nodes = x
        self.set_new_model(self.nodes)
        self.synthesize_and_compute_rf(compute_rf=False)
        self.chi2 = self.compute_chi2(only_chi2=True)
        return self.chi2

    def _callback_general(self, x):
        if (self.verbose >= 2):
            self.logger.info('chi2: {0}'.format(self.chi2))

    def invert_external(self, algorithm, use_jacobian=False, **kwargs):
        """
        Invert all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """

        for k, v in self.spectrum.items():
            v.factor_chi2 = 1.0 / (v.noise**2 * v.dof)
            
        for self.cycle in range(self.n_cycles):
            if (self.verbose >= 2):
                self.logger.info('-------------')
                self.logger.info('  Cycle {0}  '.format(self.cycle))
                
                for k, v in self.spectrum.items():
                    self.logger.info('  Weights for region {0} : SI={1} - SQ={2} - SU={3} - SV={4}'.format(k, v.stokes_weights[0,self.cycle], v.stokes_weights[1,self.cycle],
                        v.stokes_weights[2,self.cycle], v.stokes_weights[3,self.cycle]))

                self.logger.info('-------------')
                

            self.find_active_parameters(self.cycle)

            tmp = [pars['atm'] for pars in self.active_meta]
            tmp = list(set(tmp))

            self.n_free_parameters_cycle = 0
            
            for k, v in self.atmospheres.items():
                if (k in tmp):
                    if (self.verbose >= 3):
                        self.logger.info('Free parameters for {0}'.format(k))
                    for pars in self.active_meta:
                        if (pars['atm'] == k):
                            if (self.verbose >= 3):
                                if (pars['coupled'] is False):
                                    if (pars['n_nodes'] == 1):
                                        self.logger.info('  - {0} with {1} node'.format(pars['parameter'], pars['n_nodes']))
                                    else:
                                        self.logger.info('  - {0} with {1} nodes'.format(pars['parameter'], pars['n_nodes']))
                                else:
                                    self.logger.info('  - {0} coupled to {1} variable'.format(pars['parameter'], pars['n_nodes']))
                            if (pars['coupled'] is False):
                                self.n_free_parameters_cycle += pars['n_nodes']

            n_pars = len(self.nodes)

            if (use_jacobian):
                tmp = algorithm(self._func_grad, self.nodes, jac=True, callback=self._callback_general, **kwargs)
            else:
                tmp = algorithm(self._func_nograd, self.nodes, callback=self._callback_general, **kwargs)

            self._func_grad(tmp['x'])
            self._callback_general(tmp['x'])

            self.set_new_model(tmp['x'])

            self.flatten_parameters_to_reference(self.cycle)
            
    def read_observation(self):
        for k, v in self.spectrum.items():
            v.read_observation(pixel=self.pixel)
            # v.read_straylight(pixel=self.pixel)
