from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import i0_allen
from hazel.codes import hazel_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_hazel_file
import copy
# from ipdb import set_trace as stop

__all__ = ['Hazel_atmosphere']

class Hazel_atmosphere(General_atmosphere):
    def __init__(self, working_mode):
    
        super().__init__('chromosphere')

        self.height = 3.0
        self.working_mode = working_mode

        self.parameters['Bx'] = 0.0
        self.parameters['By'] = 0.0
        self.parameters['Bz'] = 0.0        
        self.parameters['tau'] = 1.0
        self.parameters['v'] = 0.0
        self.parameters['deltav'] = 8.0
        self.parameters['beta'] = 1.0
        self.parameters['a'] = 0.0
        self.parameters['ff'] = np.log(1.0)

        self.units['Bx'] = 'G'
        self.units['By'] = 'G'
        self.units['Bz'] = 'G'        
        self.units['tau'] = 'none'
        self.units['v'] = 'km/s'
        self.units['deltav'] = 'km/s'
        self.units['beta'] = 'none'
        self.units['a'] = 'none'
        self.units['ff'] = 'none'

        self.nodes['Bx'] = 0.0
        self.nodes['By'] = 0.0
        self.nodes['Bz'] = 0.0        
        self.nodes['tau'] = 0.0
        self.nodes['v'] = 0.0
        self.nodes['deltav'] = 0.0
        self.nodes['beta'] = 0.0
        self.nodes['a'] = 0.0
        self.nodes['ff'] = 0.0

        self.error['Bx'] = 0.0
        self.error['By'] = 0.0
        self.error['Bz'] = 0.0        
        self.error['tau'] = 0.0
        self.error['v'] = 0.0
        self.error['deltav'] = 0.0
        self.error['beta'] = 0.0
        self.error['a'] = 0.0
        self.error['ff'] = 0.0

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

        self.epsilon['Bx'] = 1.0
        self.epsilon['By'] = 1.0
        self.epsilon['Bz'] = 1.0      
        self.epsilon['tau'] = 1.0
        self.epsilon['v'] = 1.0
        self.epsilon['deltav'] = 1.0
        self.epsilon['beta'] = 1.0
        self.epsilon['a'] = 1.0
        self.epsilon['ff'] = 1.0

        self.jacobian['Bx'] = 1.0
        self.jacobian['By'] = 1.0
        self.jacobian['Bz'] = 1.0      
        self.jacobian['tau'] = 1.0
        self.jacobian['v'] = 1.0
        self.jacobian['deltav'] = 1.0
        self.jacobian['beta'] = 1.0
        self.jacobian['a'] = 1.0
        self.jacobian['ff'] = 1.0

        
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
        self.line_to_index = {'10830': 1, '3888': 2, '7065': 3,' 5876': 4}
        self.active_line = line
        self.wavelength_range = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top+1]

    def set_parameters(self, pars, ff):
        """
        Set the parameters of this model chromosphere

        Parameters
        ----------
        pars : list of float
            This list contains the following parameters in order: Bx, By, Bz, tau, v, delta, beta, a

        ff : float
            Filling factor

        Returns
        -------
        None
        """

        self.parameters['Bx'] = pars[0]
        self.parameters['By'] = pars[1]
        self.parameters['Bz'] = pars[2]        
        self.parameters['tau'] = pars[3]
        self.parameters['v'] = pars[4]
        self.parameters['deltav'] = pars[5]
        self.parameters['beta'] = pars[6]
        self.parameters['a'] = pars[7]
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
            

        self.model_handler = Generic_hazel_file(model_file)
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
                self.parameters[k] = self.reference[k] + np.squeeze(self.nodes[k])
            else:
                self.parameters[k] = self.reference[k]            
                            
    def print_parameters(self, first=False, error=False):
        # if (first):
        print("     {0}        {1}        {2}        {3}       {4}       {5}      {6}      {7}".format('Bx', 'By', 'Bz', 'tau', 'v', 'deltav', 'beta', 'a'), flush=True)
        print("{0:8.3f}  {1:8.3f}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.3f}  {6:8.3f}  {7:8.3f}".format(self.parameters['Bx'], self.parameters['By'], self.parameters['Bz'], self.parameters['tau'], 
            self.parameters['v'], self.parameters['deltav'], self.parameters['beta'], self.parameters['a']), flush=True)
        
        if (error):            
            print("{0:8.3f}  {1:8.3f}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.3f}  {6:8.3f}  {7:8.3f}".format(self.error['Bx'], self.error['By'], self.error['Bz'], self.error['tau'], 
            self.error['v'], self.error['deltav'], self.error['beta'], self.error['a']), flush=True)


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
                
        B = np.sqrt(self.parameters['Bx']**2 + self.parameters['By']**2 + self.parameters['Bz']**2)
        if (B == 0):
            thetaB = 0.0
        else:
            thetaB = 180.0 / np.pi * np.arccos(self.parameters['Bz'] / B)
        phiB = 180.0 / np.pi * np.arctan2(self.parameters['By'], self.parameters['Bx'])
        B1Input = np.asarray([B, thetaB, phiB])

        hInput = self.height
        tau1Input = self.parameters['tau']
        transInput = self.line_to_index[self.active_line]
        anglesInput = self.spectrum.los
        lambdaAxisInput = self.wvl_axis - self.multiplets[self.active_line]
        nLambdaInput = len(lambdaAxisInput)
                
        if (stokes is None):
            boundaryInput  = np.asfortranarray(np.zeros((4,nLambdaInput)))
            # boundaryInput[0,:] = hsra_continuum(self.multiplets[self.active_line]) #i0_allen(self.multiplets[self.active_line],1.0)            
            boundaryInput[0,:] = i0_allen(self.multiplets[self.active_line], self.spectrum.mu)
            boundaryInput *= self.spectrum.boundary
        else:            
            # boundaryInput = np.asfortranarray(stokes * hsra_continuum(self.multiplets[self.active_line])) #i0_allen(self.multiplets[self.active_line],1.0)
            boundaryInput = np.asfortranarray(stokes)
                    
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
        
        return ff * stokes #/ hsra_continuum(self.multiplets[self.active_line])