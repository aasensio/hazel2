from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import i0_allen
from hazel.codes import hazel_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_hazel_file
import copy
from ipdb import set_trace as stop

__all__ = ['Hazel_atmosphere']

class Hazel_atmosphere(General_atmosphere):
    def __init__(self):
    
        super().__init__('chromosphere')

        
    def add_active_line(self, line, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------        
        lines : str
            Line to activate: ['10830','5876']
        spectrum : Spectrum
            Spectrum object
        wvl_range : float
            Vector containing wavelength range over which to synthesize this line
        
        Returns
        -------
        None
    
        """        
        self.active_line = line
        self.wavelength_range = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top+1]

        self.height = 3.0

        self.parameters['Bx'] = 0.0
        self.parameters['By'] = 0.0
        self.parameters['Bz'] = 0.0        
        self.parameters['tau'] = 1.0
        self.parameters['v'] = 0.0
        self.parameters['deltav'] = 8.0
        self.parameters['beta'] = 1.0
        self.parameters['a'] = 0.0
        self.parameters['ff'] = np.log(1.0)

        self.nodes['Bx'] = 0.0
        self.nodes['By'] = 0.0
        self.nodes['Bz'] = 0.0        
        self.nodes['tau'] = 0.0
        self.nodes['v'] = 0.0
        self.nodes['deltav'] = 0.0
        self.nodes['beta'] = 0.0
        self.nodes['a'] = 0.0
        self.nodes['ff'] = 0.0

        self.n_nodes['Bx'] = 0
        self.n_nodes['By'] = 0
        self.n_nodes['Bz'] = 0        
        self.n_nodes['tau'] = 0
        self.n_nodes['v'] = 0
        self.n_nodes['deltav'] = 0
        self.n_nodes['beta'] = 0
        self.n_nodes['a'] = 0
        self.n_nodes['ff'] = 0

        self.ranges['Bx'] = None
        self.ranges['By'] = None
        self.ranges['Bz'] = None
        self.ranges['tau'] = None
        self.ranges['v'] = None
        self.ranges['deltav'] = None
        self.ranges['beta'] = None
        self.ranges['a'] = None
        self.ranges['ff'] = None

        self.cycles['Bx'] = None
        self.cycles['By'] = None
        self.cycles['Bz'] = None        
        self.cycles['tau'] = None
        self.cycles['v'] = None
        self.cycles['deltav'] = None
        self.cycles['beta'] = None
        self.cycles['a'] = None
        self.cycles['ff'] = None

        self.epsilon['Bx'] = 500.0
        self.epsilon['By'] = 500.0
        self.epsilon['Bz'] = 500.0      
        self.epsilon['tau'] = 1.0
        self.epsilon['v'] = 5.0
        self.epsilon['deltav'] = 5.0
        self.epsilon['beta'] = 1.0
        self.epsilon['a'] = 0.5
        self.epsilon['ff'] = 1.0

    def set_reference(self):
        """
        Set reference model to that of the current parameters

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.nodes_to_model()
        self.reference = copy.deepcopy(self.parameters)

    def set_parameters(self, pars):
        self.parameters['Bx'] = pars[0]
        self.parameters['By'] = pars[1]
        self.parameters['Bz'] = pars[2]        
        self.parameters['tau'] = pars[3]
        self.parameters['v'] = pars[4]
        self.parameters['deltav'] = pars[5]
        self.parameters['beta'] = pars[6]
        self.parameters['a'] = pars[7]
        self.parameters['ff'] = pars[8]

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
            if (verbose):
                print('    * Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            self.model_filename = model_file
            self.model_handler = Generic_hazel_file(model_file)
            self.model_handler.open()
            out = self.model_handler.read()
            self.model_handler.close()

            # out = np.loadtxt(model_file, skiprows=1)
            self.set_parameters(out)
            self.reference = copy.deepcopy(self.parameters)
        
        if (extension == 'h5'):
            if (verbose):
                print('    * Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            self.model_handler = Generic_hazel_file(model_file)
            

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
                self.parameters[k] = self.reference[k] + np.squeeze(self.nodes[k])

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

        self.nodes_to_model()
                
        B = np.sqrt(self.parameters['Bx']**2 + self.parameters['By']**2 + self.parameters['Bz']**2)
        if (B == 0):
            thetaB = 0.0
        else:
            thetaB = 180.0 / np.pi * np.arccos(self.parameters['Bz'] / B)
        phiB = 180.0 / np.pi * np.arctan2(self.parameters['By'], self.parameters['Bx'])
        B1Input = np.asarray([B, thetaB, phiB])

        hInput = self.height
        tau1Input = self.parameters['tau']
        transInput = 1
        anglesInput = np.asarray([0.0,0.0,90.0])
        lambdaAxisInput = self.wvl_axis - self.multiplets[self.active_line]
        nLambdaInput = len(lambdaAxisInput)
        
        if (stokes is None):
            boundaryInput  = np.asfortranarray(np.zeros((4,nLambdaInput)))
            boundaryInput[0,:] = hsra_continuum(self.multiplets[self.active_line]) #i0_allen(self.multiplets[self.active_line],1.0)            
        else:            
            boundaryInput = np.asfortranarray(stokes * hsra_continuum(self.multiplets[self.active_line])) #i0_allen(self.multiplets[self.active_line],1.0)
                    
        dopplerWidthInput = self.parameters['deltav']
        dampingInput = self.parameters['a']
        dopplerVelocityInput = self.parameters['v']
        betaInput = self.parameters['beta']
        nbarInput = np.asarray([0.0,0.0,0.0,0.0])
        omegaInput = np.asarray([0.0,0.0,0.0,0.0])
        
        args = (self.index, B1Input, hInput, tau1Input, boundaryInput, transInput, 
            anglesInput, nLambdaInput, lambdaAxisInput, dopplerWidthInput, 
            dampingInput, dopplerVelocityInput, 
            betaInput, nbarInput, omegaInput)
        
        l, stokes = hazel_code._synth(*args)

        ff = self.parameters['ff']
        
        return ff * stokes / hsra_continuum(self.multiplets[self.active_line])