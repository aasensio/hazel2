from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import fvoigt, i0_allen
from hazel.hsra import hsra_continuum
from hazel.io import Generic_stray_file
import copy
import scipy.constants as constants

__all__ = ['Straylight_atmosphere']

class Straylight_atmosphere(General_atmosphere):
    def __init__(self, working_mode, name=''):
    
        super().__init__('straylight', name=name)

        self.working_mode = working_mode

        self.parameters['ff'] = 1.0
        self.parameters['v'] = 0.0

        self.nodes['ff'] = 0.0
        self.nodes['v'] = 0.0
        
        self.n_nodes['ff'] = 0
        self.n_nodes['v'] = 0
        
        self.epsilon['ff'] = 1.0
        self.epsilon['v'] = 1.0
        
        self.ranges['ff'] = None
        self.ranges['v'] = None
        
        self.cycles['ff'] = None
        self.cycles['v'] = None
        

    def add_active_line(self, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------        
        None

        Returns
        -------
        None
        """

        self.wavelength_range = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top+1]
        self.avg_wavelength = np.mean(self.wvl_axis)
        self.stokes = np.zeros((4,len(self.wvl_axis)))
        
        
    def set_parameters(self, pars, ff):
        self.stray_profile = pars[0]
        self.parameters['v'] = pars[1]
        self.parameters['ff'] = ff

        # Check that parameters are inside borders by clipping inside the interval with a border of 1e-8
        if (self.working_mode == 'inversion'):
            for k, v in self.parameters.items():            
                self.parameters[k] = np.clip(v, self.ranges[k][0] + 1e-8, self.ranges[k][1] - 1e-8)
    
    def set_straylight(self, stokes):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            verbosity flag

        Returns
        -------
            None
        """
        pass
        # self.stray_stokesI = self.spectrum.

    def load_reference_model(self, model_file, verbose):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            verbosity flag

        Returns
        -------
            None
        """
        extension = os.path.splitext(model_file)[1][1:]
        if (extension == '1d'):
            if (verbose >= 1):
                self.logger.info('    * Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            self.model_filename = model_file
                    
        if (extension == 'h5'):
            if (verbose >= 1):
                self.logger.info('    * Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'            

        self.model_handler = Generic_stray_file(model_file)
        self.model_handler.open()
        out, ff = self.model_handler.read(pixel=0)
        self.model_handler.close()
        
        self.set_parameters(out, ff)

        self.init_reference()

    def nodes_to_model(self):
        """
        Transform from nodes to model
        
        Parameters
        ----------
        None
                                
        Returns
        -------
        None
        """        
        for k, v in self.nodes.items():
            if (self.n_nodes[k] > 0):
                self.parameters[k] = self.reference[k] + self.nodes[k]
            else:
                self.parameters[k] = self.reference[k]

    def print_parameters(self, first=False, error=False):
        self.logger.info("     {0}        {1}".format('v', 'ff'))

        self.logger.info("     {0:8.3f}  {1:8.3f}  ".format(np.atleast_1d(self.parameters['v'])[0], \
            np.atleast_1d(self.parameters['ff'])[0]))

    def synthesize(self, nlte=None):
        """
        Carry out the synthesis and returns the Stokes parameters
        
        Parameters
        ----------
        None
                        
        Returns
        -------        
        stokes : float
            Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (4,nLambda)        
        """
        
        if (self.working_mode == 'inversion'):
            self.nodes_to_model()
            self.to_physical()
        
        dlambda = self.parameters['v'] * 1e5 / (100.0 * constants.c) * self.avg_wavelength
        
        # self.stokes[0,:] = np.interp(self.wvl_axis - dlambda, self.wvl_axis, self.stray_profile[self.wvl_range[0]:self.wvl_range[1]])
        self.stokes[0,:] = np.interp(self.wvl_axis - dlambda, self.wvl_axis, self.stray_profile)
    
        error = 0
                                        
        return self.parameters['ff'] * self.stokes * i0_allen(np.mean(self.wvl_axis), self.spectrum.mu), error #* hsra_continuum(np.mean(self.wvl_axis)), error