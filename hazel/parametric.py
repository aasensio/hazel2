from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import fvoigt
from hazel.hsra import hsra_continuum
from hazel.io import Generic_parametric_file
import copy
# from ipdb import set_trace as stop

__all__ = ['Parametric_atmosphere']

class Parametric_atmosphere(General_atmosphere):
    def __init__(self, working_mode):
    
        super().__init__('parametric')

        self.working_mode = working_mode

        self.parameters['lambda0'] = 0.0
        self.parameters['sigma'] = 0.0
        self.parameters['depth'] = 0.0
        self.parameters['a'] = 0.0
        self.parameters['ff'] = np.log(1.0)

        self.nodes['lambda0'] = 0.0
        self.nodes['sigma'] = 0.0
        self.nodes['depth'] = 0.0
        self.nodes['a'] = 0.0
        self.nodes['ff'] = 0.0

        self.n_nodes['lambda0'] = 0
        self.n_nodes['sigma'] = 0
        self.n_nodes['depth'] = 0
        self.n_nodes['a'] = 0
        self.n_nodes['ff'] = 0

        self.epsilon['lambda0'] = 1.0
        self.epsilon['sigma'] = 0.1
        self.epsilon['depth'] = 0.5
        self.epsilon['a'] = 0.05
        self.epsilon['ff'] = 0.05

        self.ranges['lambda0'] = None
        self.ranges['sigma'] = None
        self.ranges['depth'] = None
        self.ranges['a'] = None
        self.ranges['ff'] = None
        
        self.cycles['lambda0'] = None
        self.cycles['sigma'] = None
        self.cycles['depth'] = None        
        self.cycles['a'] = None
        self.cycles['ff'] = None
        

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
    
    def set_parameters(self, pars, ff):
        """
        Set parameters of the model

        Parameters
        ----------
        pars : list of floats
            Values of the parameters: lambda0, sigma, depth, a

        ff : float
            Filling factor

        Returns
        -------
        None
        """
        self.parameters['lambda0'] = pars[0]
        self.parameters['sigma'] = pars[1]
        self.parameters['depth'] = pars[2]
        self.parameters['a'] = pars[3]
        self.parameters['ff'] = ff

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
                print('    * Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            self.model_filename = model_file
            
        
        if (extension == 'h5'):
            if (verbose >= 1):
                print('    * Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            
        self.model_handler = Generic_parametric_file(model_file)
        self.model_handler.open()
        out, ff = self.model_handler.read()
        self.model_handler.close()

        # out = np.loadtxt(model_file, skiprows=1)
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

    
    def synthesize(self, stokes=None):
        """
        Carry out the synthesis and returns the Stokes parameters
        
        Parameters
        ----------
        stokes : float
            An array of size [4 x nLambda] with the input Stokes parameter.
                        
        Returns
        -------        
        stokes : float
            Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (4,nLambda)        
        """
        
        if (self.working_mode == 'inversion'):
            self.nodes_to_model()
            self.to_physical()

        lambda0 = self.parameters['lambda0']
        sigma = self.parameters['sigma']
        d = self.parameters['depth']
        a = self.parameters['a']
        ff = self.parameters['ff']
                
        if (stokes is None):
            n_lambda = len(self.wvl_axis)
            stokes = np.zeros((4,n_lambda))
            stokes[0,:] = hsra_continuum(lambda0)

        profile, _ = fvoigt(a, (self.wvl_axis - lambda0) / sigma)

        stokes[0,:] *= ff * (1.0 - d * profile)
        
        return stokes