from collections import OrderedDict
import numpy as np
import os
from hazel.util import i0_allen, fvoigt
from hazel.codes import hazel_code, sir_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_hazel_file, Generic_SIR_file, Generic_parametric_file
from hazel.transforms import transformed_to_physical, physical_to_transformed, jacobian_transformed_to_physical, jacobian_transformation
import copy
from hazel.sir import Sir
import logging

# from ipdb import set_trace as stop


__all__ = ['General_atmosphere']
    
class General_atmosphere(object):
    def __init__(self, atm_type, name):
        self.ff = 1.0
        self.name = name

        self.logger = logging.getLogger("model")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []
        ch = logging.StreamHandler()        
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        self.active_lines = []
        self.wavelength = dict()
        self.wavelength_range = dict()
        self.wvl_axis = dict()
        self.wvl_range = dict()
        self.type = atm_type
        self.spectrum = dict()
        self.active = False
        self.n_pixel = 1
        self.stray_profile = None
        

        self.multiplets = {'10830': 10829.0911, '3888': 3888.6046, '7065': 7065.7085, '5876': 5875.9663}

        self.parameters = OrderedDict()
        self.nodes_location = OrderedDict()
        self.ranges = OrderedDict()
        self.regularization = OrderedDict()
        self.cycles = OrderedDict()
        self.n_nodes = OrderedDict()
        self.nodes = OrderedDict()
        self.epsilon = OrderedDict()
        self.error = OrderedDict()
        self.jacobian = OrderedDict()
        self.units = OrderedDict()

        self.eps_borders = 1e-5

    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d:
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d:
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)

    def __str__(self):
        tmp = ''
        for l, par in self.parameters.items():
            tmp += '{0}: {1} {2}\n'.format(l, par, self.units[l])
        return tmp

    def allocate_info_cycles(self, n_cycles):
        """
        Set the appropriate variables to store per-cycle models
        
        Parameters
        ----------        
        n_cycles : int
            Number of cycles
        
        Returns
        -------
        None
    
        """

        self.reference_cycle = [None] * n_cycles
        self.error_cycle = [None] * n_cycles
        self.nodes_location_cycle = [None] * n_cycles

    def to_physical(self):
        """
        Transform the atmospheric parameters from transformed domain to physical domain given the ranges.
        This only applies in inversion mode
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """
        
        for k, v in self.parameters.items():            
            if (k in self.ranges):
                lower = self.ranges[k][0] #- self.eps_borders
                upper = self.ranges[k][1] #+ self.eps_borders
                self.parameters[k] = transformed_to_physical(v, lower, upper)
                        
    def to_transformed(self):
        """
        Transform the atmospheric parameters from transformed domain to physical domain given the ranges.
        This only applies in inversion mode
        
        Parameters
        ----------        
        None
        
        Returns
        -------
        None
    
        """
        for k, v in self.parameters.items():
            if (k in self.ranges):
                lower = self.ranges[k][0] #- self.eps_borders
                upper = self.ranges[k][1] #+ self.eps_borders            
                self.parameters[k] = physical_to_transformed(v, lower, upper)
            

    def set_reference(self, cycle=None):
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
        self.nodes_location =  copy.deepcopy(self.nodes_location)

        if (cycle is not None):
            self.to_physical()        
            self.reference_cycle[cycle] = copy.deepcopy(self.parameters)
            self.nodes_location_cycle[cycle] =  copy.deepcopy(self.nodes_location)
            self.error_cycle[cycle] = copy.deepcopy(self.error)

    def reset_reference(self):
        """
        Reset reference model to the original reference loaded from file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """                    
        self.reference = copy.deepcopy(self.original_reference)

    def init_reference(self, check_borders=False):
        """
        Initialize the reference atmosphere to the values of the parameters, doing the inverse transformation if in inversion mode

        Parameters
        ----------
        check_borders : bool
            Check that the input parameters are inside the ranges of parameters

        Returns
        -------
        None
        """

        # Transform parameters to  unbounded domain
        if (self.working_mode == 'inversion'):
            if (check_borders):
                for k, v in self.parameters.items():                    
                    if (k in self.ranges):
                        if (not np.all(np.logical_and(v >= self.ranges[k][0], v <= self.ranges[k][1]))):
                            raise Exception("Parameter {0} of atmosphere {1} is outside ranges".format(k, self.name))

            self.to_transformed()            

        self.reference = copy.deepcopy(self.parameters)
        self.original_reference = copy.deepcopy(self.parameters)        