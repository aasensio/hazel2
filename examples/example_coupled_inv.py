import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py

# Test a single inversion in non-iterator mode without randomization
mod = hazel.Model('configurations/conf_inv_D3_10830.ini', working_mode='inversion', verbose=4)
mod.read_observation()
mod.open_output()

mod.invert()
mod.write_output()

mod.close_output()

fig, ax = pl.subplots(nrows=1, ncols=2)

f = h5py.File('output.h5', 'r')
obs1 = np.loadtxt('observations/multi_10830_stokes.1d', skiprows=7)
obs2 = np.loadtxt('observations/multi_D3_stokes.1d', skiprows=7)

ax[0].plot(obs1[:,0])
ax[0].plot(f['spec1']['stokes'][0,0,0,0,:])

ax[1].plot(obs2[:,0])
ax[1].plot(f['spec2']['stokes'][0,0,0,0,:])

pl.show()

f.close()
