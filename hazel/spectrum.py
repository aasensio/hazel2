import numpy as np
from hazel.io import Generic_observed_file, Generic_stray_file
# from ipdb import set_trace as stop

__all__ = ['Spectrum']

class Spectrum(object):
    def __init__(self, wvl=None, weights=None, observed_file=None, stray=None, name=None, stokes_weights=None):
        
        self.wavelength_axis = None
        self.stokes = None
        self.stokes_perturbed = None
        self.pixel = 0

        if (wvl is not None):
            self.add_spectrum(wvl)

        if (weights is not None):
            self.add_weights(weights)

        if (observed_file is not None):
            self.add_observed_file(observed_file)

        if (stray is not None):
            self.add_stray_file(stray)
            self.read_straylight()

        if (name is not None):
            self.add_name(name)

        if (stokes_weights is not None):
            self.add_stokes_weights(stokes_weights)

    def add_spectrum(self, wvl):
        """
        Add a new spectrum, initializing the Stokes parameters containers
        
        Parameters
        ----------        
        wvl : float
            Wavelength axis
        
        Returns
        -------
        None
    
        """  
        self.wavelength_axis = wvl
        self.stokes = np.zeros((4,len(wvl)))
        self.stokes_perturbed = np.zeros((4,len(wvl)))
        self.stray = np.zeros((4,len(wvl)))
        self.obs = np.zeros((4,len(wvl)))
        self.noise = np.zeros((4,len(wvl)))
        self.dof = 4.0 * len(wvl)

    def add_weights(self, weights):
        """
        Add new weights for the Stokes parameters
        
        Parameters
        ----------        
        weights : float
            Weights
        
        Returns
        -------
        None
    
        """  
        self.wavelength_weights = weights

    def add_observed_file(self, observed_file):        
        """
        Add a new file with observations and open it
        
        Parameters
        ----------        
        observed_file : str
            File with the observations
        
        Returns
        -------
        None
    
        """  
        self.observed_handle = Generic_observed_file(observed_file)
        self.n_pixel = self.observed_handle.get_npixel()

    def add_stray_file(self, stray_file):        
        """
        Add a new file with straylight and open it
        
        Parameters
        ----------        
        straylight_file : str
            File with the straylight
        
        Returns
        -------
        None
    
        """  
        self.stray_handle = Generic_stray_file(stray_file)
        self.stray_handle.open()

    def add_name(self, name):
        """
        Add a new name to this spectral region
        
        Parameters
        ----------        
        name : str
            Name of the region
        
        Returns
        -------
        None
    
        """  
        self.name = name

    def add_stokes_weights(self, weights):
        """
        Add Stokes weights
        
        Parameters
        ----------        
        weights : float
            Array of size 4 with the weights for all Stokes parameters
        
        Returns
        -------
        None
    
        """  
        self.stokes_weights = weights

    def next_pixel(self):
        """
        Skip to next pixel
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """                  
        self.pixel += 1

    def open_observation(self):
        """
        Open the next pixel with observations
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """          
        self.observed_handle.open()

    def read_observation(self, pixel=None):
        """
        Read the next pixel with observations
        
        Parameters
        ----------        
        pixel : int
            Pixel to read
        
        Returns
        -------
        None
    
        """          
        self.obs, self.noise = self.observed_handle.read(pixel=pixel)
        stop()

    def read_straylight(self, pixel=None):
        """
        Read the next pixel with straylight
        
        Parameters
        ----------        
        pixel : int
            Pixel to read
        
        Returns
        -------
        None
    
        """          
        self.stray = self.stray_handle.read(pixel=pixel)