import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py

# Test a single inversion in non-iterator mode without randomization
mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion', verbose=2)
mod.read_observation()
mod.open_output()

mod.invert()
mod.write_output()

mod.close_output()

final = np.loadtxt('photospheres/model_photosphere.1d', skiprows=4)
start = np.loadtxt('photospheres/model_photosphere_200.1d', skiprows=4)
f = h5py.File('output.h5')
fig = pl.figure()
pl.plot(f['ph1']['T'][0,0,0,:])
pl.plot(final[:,1])
pl.plot(start[:,1], 'x')
pl.show()
f.close()

# Test a single inversion in non-iterator mode
# mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion', verbose=2, randomization=2)
# mod.read_observation()
# mod.open_output()

# for loop in range(2):
#     mod.invert(randomize=True)
#     mod.write_output(randomization=loop)

# mod.close_output()

# final = np.loadtxt('photospheres/model_photosphere.1d', skiprows=4)
# start = np.loadtxt('photospheres/model_photosphere_200.1d', skiprows=4)
# f = h5py.File('output.h5')
# fig = pl.figure()
# pl.plot(f['ph1']['T'][0,0,0,:])
# pl.plot(f['ph1']['T'][0,1,0,:])
# pl.plot(final[:,1])
# pl.plot(start[:,1], 'x')
# pl.show()
# f.close()