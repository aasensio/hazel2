import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from scipy.optimize import minimize
import gc

# Test a single inversion in non-iterator mode
mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion', verbose=2)
mod.read_observation()
mod.open_output()
mod.invert()
mod.write_output()
mod.close_output()


final = np.loadtxt('photospheres/model_photosphere.1d', skiprows=4)
start = np.loadtxt('photospheres/model_photosphere_200.1d', skiprows=4)
f = h5py.File('output.h5', 'r')
pl.plot(f['ph1']['T'][0,0,:], label='inverted')
pl.plot(final[:,1], label='target')
pl.plot(start[:,1], 'x', label='initial')
f.close()
pl.legend()

mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion', verbose=2)
mod.read_observation()
mod.open_output()
mod.invert_external(minimize, method='Nelder-Mead')
mod.write_output()
mod.close_output()

final = np.loadtxt('photospheres/model_photosphere.1d', skiprows=4)
start = np.loadtxt('photospheres/model_photosphere_200.1d', skiprows=4)
f = h5py.File('output.h5', 'r')

pl.figure()
pl.plot(f['ph1']['T'][0,0,:], label='inverted')
pl.plot(final[:,1], label='target')
pl.plot(start[:,1], 'x', label='initial')
f.close()
pl.legend()

pl.show()
