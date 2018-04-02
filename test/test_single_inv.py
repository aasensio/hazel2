import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

# Test a single inversion in non-iterator mode
mod = hazel.Model('conf_single.ini', working_mode='inversion')
mod.read_observation()
mod.open_output()
mod.invert()
mod.close_output()