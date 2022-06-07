import numpy as np
from hazel.io import Generic_observed_file, Generic_stray_file, Generic_mask_file
import os.path
from astropy.constants import c
from scipy import interpolate

#EDGAR: remember here the routine add_boundary

__all__ = ['Spectrum']

class Spectrum(object):
    def __init__(self, wvl=None, weights=None, observed_file=None, name=None, stokes_weights=None, los=None, boundary=None, mask_file=None, 
        instrumental_profile=None, save_all_cycles=False, root='', wvl_lr=None):
        
        self.wavelength_axis = None
        self.stokes = None
        self.stokes_perturbed = None
        self.pixel = 0
        self.boundary_single = boundary
        self.psf_spectral = None
        self.chi2 = -1.0
        self.rss = -1.0
        self.normalization = 'on-disk'
        self.save_all_cycles = save_all_cycles
        self.root = root
        
        if (wvl is not None):
            self.add_spectrum(wvl, wvl_lr)

        if (weights is not None):
            self.add_weights(weights)

        if (observed_file is not None):
            self.add_observed_file(observed_file)

        self.add_mask_file(mask_file)

        # if (stray is not None):
            # self.add_stray_file(stray)
            # self.read_straylight()

        if (name is not None):
            self.add_name(name)

        if (stokes_weights is not None):
            self.add_stokes_weights(stokes_weights)

        if (los is not None):
            self.los = los
            self.mu = np.cos(self.los[0] * np.pi / 180.0)

        if (instrumental_profile is not None):

            n = len(self.wavelength_axis)

            # Instrumental profile given as a file
            if os.path.exists(self.root + instrumental_profile):
                tmp = np.loadtxt(self.root + instrumental_profile)
                f = interpolate.interp1d(tmp[:,0] + self.wavelength_axis[n//2-1], tmp[:,1], kind='cubic', fill_value=0.0, bounds_error=False)                
                self.psf_spectral = f(self.wavelength_axis)
            else:
                sigma = float(instrumental_profile) * np.mean(self.wavelength_axis) / c.to('km/s').value                
                self.psf_spectral = np.exp(-(self.wavelength_axis - self.wavelength_axis[n//2-1])**2 / sigma**2)

            self.psf_spectral /= np.sum(self.psf_spectral)
                
    def allocate_info_cycles(self, n_cycles):
        """
        Set the appropriate variables to store per-cycle spectra
        
        Parameters
        ----------        
        n_cycles : int
            Number of cycles
        
        Returns
        -------
        None
    
        """

        self.stokes_cycle = [None] * n_cycles
        self.stokes_lr_cycle = [None] * n_cycles
        self.chi2_cycle = [-1.0] * n_cycles
        self.bic_cycle = [-1.0] * n_cycles
        self.aic_cycle = [-1.0] * n_cycles

    def add_spectrum(self, wvl, wvl_lr):
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
        self.wavelength_axis_lr = wvl_lr
        self.stokes = np.zeros((4,len(wvl)))
        self.stokes_perturbed = np.zeros((4,len(wvl)))
        self.stray = np.zeros((4,len(wvl)))
        self.obs = np.zeros((4,len(wvl)))
        self.noise = np.zeros((4,len(wvl)))
        self.dof = 4.0 * len(wvl)
        if (self.boundary_single is not None):
            
            self.boundary = self.boundary_single[:,None] * np.ones((1,len(wvl)))

            if (self.boundary_single[0] == 0.0):                
                self.normalization = 'off-limb'
            else:
                self.normalization = 'on-disk'

        if (self.wavelength_axis_lr is not None):
            self.interpolate_to_lr = True
            self.stokes_lr = np.zeros((4,len(wvl_lr)))
            self.stokes_perturbed_lr = np.zeros((4,len(wvl_lr)))
        else:
            self.interpolate_to_lr = False
            
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
        self.observed_handle = Generic_observed_file(observed_file, self.root)
        self.n_pixel = self.observed_handle.get_npixel()

    def add_mask_file(self, mask_file):        
        """
        Add a new file with a mask
        
        Parameters
        ----------        
        mask_file : str
            File with the mask
        
        Returns
        -------
        None
    
        """  
        if (mask_file is not None):
            self.mask_handle = Generic_mask_file(mask_file, self.root)            
            n_pixel = self.mask_handle.get_npixel()
            if (n_pixel != self.n_pixel):
                raise Exception("Number of pixels in mask is different from number of observed pixels")
        else:
            self.mask_handle = Generic_mask_file(None, None)
        
    # def add_stray_file(self, stray_file):        
    #     """
    #     Add a new file with straylight and open it
        
    #     Parameters
    #     ----------        
    #     straylight_file : str
    #         File with the straylight
        
    #     Returns
    #     -------
    #     None
    
    #     """
    #     self.stray_present = True
    #     self.stray_handle = Generic_stray_file(stray_file)
    #     self.stray_handle.open()

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

    def set_los(self, los):
        """
        Set a new value for the LOS
        
        Parameters
        ----------        
        los : list or array of size 3
            theta, phi and gamma of the line-of-sight
        
        Returns
        -------
        None
    
        """          
        self.los = los
        self.mu = np.cos(self.los[0] * np.pi / 180.0)

    def set_boundary(self, boundary):
        """
        Set a new value for the boundary condition
        
        Parameters
        ----------        
        los : list or array of size 4
            I0, Q0, U0, V0
        
        Returns
        -------
        None
    
        """          
        self.boundary = boundary[:,None] * np.ones((1,len(self.wavelength_axis)))
        if (np.all(self.boundary[0,:]) == 0.0):
            self.normalization = 'off-limb'
        else:
            self.normalization = 'on-disk'

    def set_normalization(self, normalization):
        """
        Set a new value for the LOS
        
        Parameters
        ----------        
        los : str
            'off-limb' or 'on-disk'
        
        Returns
        -------
        None
    
        """          
        self.normalization = normalization

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
        Open the file with the observations for this spectrum
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """          
        self.observed_handle.open()

    def close_observation(self):
        """
        Close the file with the observations for this spectrum
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """          
        self.observed_handle.close()

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
        self.obs, self.noise, self.los, self.mu, self.boundary = self.observed_handle.read(pixel=pixel)
        if (np.all(self.boundary[0,:]) == 0.0):
            self.normalization = 'off-limb'
        else:
            self.normalization = 'on-disk'

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

    def read_mask(self, pixel=None):
        """
        Read the next pixel with mask
        
        Parameters
        ----------        
        pixel : int
            Pixel to read
        
        Returns
        -------
        None
    
        """                  
        if (self.mask_handle is not None):
            self.mask = self.mask_handle.read(pixel=pixel)
        else:
            self.mask = 1