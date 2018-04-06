from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import i0_allen
from hazel.codes import sir_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_SIR_file
import scipy.interpolate as interp
# from ipdb import set_trace as stop
import copy

__all__ = ['SIR_atmosphere']

# sir_parameters = OrderedDict.fromkeys('T B thetaB phiB v')

class SIR_atmosphere(General_atmosphere):
    def __init__(self):
        
        super().__init__('photosphere')
        self.ff = 1.0        
        self.macroturbulence = np.zeros(1)
        
        self.parameters['T'] = None
        self.parameters['vmic'] = None
        self.parameters['v'] = None
        self.parameters['Bx'] = None
        self.parameters['By'] = None
        self.parameters['Bz'] = None
        self.parameters['ff'] = None

        self.n_nodes['T'] = 0
        self.n_nodes['vmic'] = 0
        self.n_nodes['v'] = 0
        self.n_nodes['Bx'] = 0
        self.n_nodes['By'] = 0
        self.n_nodes['Bz'] = 0
        self.n_nodes['ff'] = 0

        self.nodes['T'] = 0
        self.nodes['vmic'] = 0
        self.nodes['v'] = 0
        self.nodes['Bx'] = 0
        self.nodes['By'] = 0
        self.nodes['Bz'] = 0
        self.nodes['ff'] = 0

        self.ranges['T'] = None
        self.ranges['vmic'] = None
        self.ranges['v'] = None
        self.ranges['Bx'] = None
        self.ranges['By'] = None
        self.ranges['Bz'] = None
        self.ranges['ff'] = None

        self.cycles['T'] = None
        self.cycles['vmic'] = None
        self.cycles['v'] = None
        self.cycles['Bx'] = None
        self.cycles['By'] = None
        self.cycles['Bz'] = None
        self.cycles['ff'] = None

        self.epsilon['T'] = 1000.0
        self.epsilon['vmic'] = 1.0
        self.epsilon['v'] = 1.0
        self.epsilon['Bx'] = 200.0
        self.epsilon['By'] = 200.0
        self.epsilon['Bz'] = 200.0
        self.epsilon['ff'] = 1.0
        
    def list_lines(self):
        """
        List the lines available in SIR for synthesis
            
        """
        f = open('LINEAS', 'r')
        lines = f.readlines()
        f.close()

        print("Available lines:")
        for l in lines[:-1]:
            print(l[:-1])

    def add_active_line(self, lines, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------
        lines : str
            Line to activate
        spectrum : Spectrum
            Spectrum object
        wvl_range : float
            Vector containing wavelength range over which to synthesize this line
        
        Returns
        -------
        None
    
        """
        
        self.lines = lines
        self.wvl_range_lambda = wvl_range

        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = np.array([ind_low, ind_top+1])
        
    def interpolate_nodes(self, log_tau, reference, nodes):
        """
        Generate a model atmosphere by interpolating the defined nodes. The interpolation
        order depends on the number of nodes.
        
        Parameters
        ----------
        log_tau : float
            Vector of log optical depth at 500 nm
        reference : float
            Vector with the reference atmosphere to which the nodes are added
        nodes : float
            List with the position of the nodes

        Returns
        -------
        real
            Vector with the interpolated atmosphere
    
        """
        n_nodes = len(nodes)
        n_depth = len(log_tau)

        if (n_nodes == 0):
            return reference
        
        if (n_nodes == 1):
            return reference + nodes[0]
        
        if (n_nodes == 2):
            pos = np.linspace(0, n_depth-1, n_nodes, dtype=int)
            # coeff = np.polyfit(log_tau[pos], nodes, 1)            
            f = interp.interp1d(log_tau[pos], nodes, 'linear')
            return reference + f(log_tau) #np.polyval(coeff, log_tau)

        if (n_nodes == 3):
            pos = np.linspace(0, n_depth-1, n_nodes, dtype=int)
            # coeff = np.polyfit(log_tau[pos], nodes, 2)      
            f = interp.interp1d(log_tau[pos], nodes, 'quadratic')
            return reference + f(log_tau)

        if (n_nodes > 3):            
            pos = np.linspace(n_depth-1, 0, n_nodes, dtype=int)
            f = interp.PchipInterpolator(log_tau[pos], nodes, extrapolate=True)            
            return reference + f(log_tau)

    def load_reference_model(self, model_file, verbose):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            Verbosity

        Returns
        -------
        None
        """
        extension = os.path.splitext(model_file)[1][1:]
        if (extension == '1d'):
            if (verbose):
                print('    * Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            self.model_handler = Generic_SIR_file(model_file)
            self.model_handler.open()
            out, ff = self.model_handler.read()
            self.model_handler.close()
            
            if (np.min(out[:,2]) > 0.0):
                self.pe_present = True
            else:
                self.pe_present = False
                        
            self.set_parameters([out, ff])
            self.reference = copy.deepcopy(self.parameters)
        
        if (extension == 'h5'):
            if (verbose):
                print('    * Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            self.model_handler = Generic_SIR_file(model_file)
            self.pe_present = True
            # self.model_file = model_file

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
            
    def set_parameters(self, model_in):
        """
        Set the parameters of the current model to those passed as argument

        Parameters
        ----------
        model_in : float
            Array with the model        
        ff : float
            Value of the filling factor

        Returns
        -------
        None
        """

        model = model_in[0]
        ff = model_in[1]

        self.parameters['log_tau'] = model[:,0]
        self.parameters['T'] = model[:,1]
        self.parameters['vmic'] = model[:,3]
        self.parameters['v'] = model[:,4]
        self.parameters['Bx'] = model[:,5]
        self.parameters['By'] = model[:,6]
        self.parameters['Bz'] = model[:,7]

        if (self.pe_present):
            self.parameters['Pe'] = model[:,2]
        else:
            self.parameters['Pe'] = -np.ones(len(self.parameters['log_tau']))
            self.parameters['Pe'][-1] = 1.11634e-1        
        
        self.parameters['ff'] = ff

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
            if (k is not 'log_tau'):
                if (self.n_nodes[k] > 0):
                    self.parameters[k] = self.interpolate_nodes(self.parameters['log_tau'], self.reference[k], self.nodes[k])
                    
    def model_to_nodes(self):
        """
        Transform from model to nodes
        
        Parameters
        ----------
        None
                                
        Returns
        -------
        None        
        """

        pass

        # for k, v in self.parameters.items():
            # if (k is not 'log_tau'):

        # self.interpolate_nodes(self.parameters['log_tau'], self.reference['T'], [1000])
        # stop()
        # for k, v in self.parameters.items():
            # if (k is not 'log_tau'):
                
                

    def synthesize(self, stokes_in, returnRF=False):
        """
        Carry out the synthesis and returns the Stokes parameters and the response 
        functions to all physical variables at all depths
        
        Parameters
        ----------
        stokes_in : float
            An array of size [4 x nLambda] with the input Stokes parameter. It is irrelevant in this case
            because we assume that all SIR atmospheres have the Planck function as boundary.
                
        returnRF : bool, optional
            Return response functions
        
        Returns
        -------
        
        stokes : float
            Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (5,nLambda)
        rf: float (optional)
            Response functions to T, Pe, vmic, B, v, theta, phi, all of size (4,nLambda,nDepth), plus the RF to macroturbulence of size (4,nLambda)
                            It is not returned if returnRF=False
        """

        self.nodes_to_model()
        
        if (returnRF):
            stokes, rf = sir_code.synthRF(self.index, self.n_lambda, self.parameters['log_tau'], self.parameters['T'], 
                self.parameters['Pe'], self.parameters['vmic'], self.parameters['v'], self.parameters['Bx'], self.parameters['By'], 
                self.parameters['Bz'], self.macroturbulence)
            return self.parameters['ff'] * stokes[1:,:], rf
        else:
            stokes = sir_code.synth(self.index, self.n_lambda, self.parameters['log_tau'], self.parameters['T'], 
                self.parameters['Pe'], self.parameters['vmic'], self.parameters['v'], self.parameters['Bx'], self.parameters['By'], 
                self.parameters['Bz'], self.macroturbulence[0])
            
            return self.parameters['ff'] * stokes[1:,:] * hsra_continuum(np.mean(self.wvl_axis))