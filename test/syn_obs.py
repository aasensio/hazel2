import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

mod = hazel.Model('conf_synobs.ini')
mod.synthesize()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()
for i in range(4):
    ax[i].plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[i,:])
pl.show()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()
for i in range(4):
    ax[i].plot(mod.spectrum['spec2'].wavelength_axis, mod.spectrum['spec2'].stokes[i,:])
pl.show()

pl.pause(0.001)

stray = np.zeros(150)
stray[20:50] = 1.0

noise = 1e-4 * np.ones((150,4))
spec = mod.spectrum['spec1'].stokes.T + 1e-4 * np.random.randn(150,4)
np.savetxt('observations/10830.wavelength', mod.spectrum['spec1'].wavelength_axis, header='lambda')
np.savetxt('observations/10830_stokes.1d', np.hstack([spec, noise]), header='SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV')
np.savetxt('observations/10830_stray.1d', np.hstack([stray]), header='SI')
np.savetxt('observations/10830.weights', np.array([1.0,1.0,1.0,1.0]), header='wI wQ wU wV')

noise = 1e-4 * np.ones((150,4))
spec = mod.spectrum['spec2'].stokes.T + 1e-4 * np.random.randn(150,4)
np.savetxt('observations/6302.wavelength', mod.spectrum['spec2'].wavelength_axis, header='lambda')
np.savetxt('observations/6302_stokes.1d', np.hstack([spec, noise]), header='SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV')
np.savetxt('observations/6302.weights', np.array([1.0,1.0,1.0,1.0]), header='wI wQ wU wV')
stop()