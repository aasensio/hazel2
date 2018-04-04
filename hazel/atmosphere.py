from collections import OrderedDict
import numpy as np
import os
from hazel.util import i0_allen, fvoigt
from hazel.codes import hazel_code, sir_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_hazel_file, Generic_SIR_file, Generic_parametric_file
# from ipdb import set_trace as stop
from hazel.sir import Sir


__all__ = ['General_atmosphere']
    
class General_atmosphere(object):
    def __init__(self, atm_type):
        self.ff = 1.0

        self.active_lines = []
        self.wavelength = dict()
        self.wavelength_range = dict()
        self.wvl_axis = dict()
        self.wvl_range = dict()
        self.type = atm_type
        self.spectrum = dict()
        self.active = False
        self.n_pixel = 1

        self.multiplets = {'10830': 10829.0911, '3888': 3888.6046, '7065': 7065.7085, '5876': 5875.9663}

        self.parameters = OrderedDict()
        self.ranges = OrderedDict()
        self.cycles = OrderedDict()
        self.n_nodes = OrderedDict()
        self.nodes = OrderedDict()
        self.epsilon = OrderedDict()