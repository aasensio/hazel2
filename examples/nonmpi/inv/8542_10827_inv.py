import numpy as np
import matplotlib.pyplot as pl
import hazel

tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1, default='hsra')
tmp.save('hsra')


mod = hazel.Model('conf_8542_10827_inv.ini', working_mode='inversion', verbose=3)
mod.read_observation()
mod.open_output()
mod.invert()
mod.write_output()
mod.close_output()

obs_ca = np.loadtxt('8542.1d', skiprows=7)
wvl_ca = np.loadtxt('8542.wavelength', skiprows=1)
obs_si = np.loadtxt('10827.1d', skiprows=7)
wvl_si = np.loadtxt('10827.wavelength', skiprows=1)

fig, ax = pl.subplots(nrows=1, ncols=2, figsize=(12, 6))

ax[0].plot(wvl_ca, obs_ca[:, 0])
ax[0].plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[0,:])

ax[1].plot(wvl_si, obs_si[:, 0])
ax[1].plot(mod.spectrum['spec2'].wavelength_axis, mod.spectrum['spec2'].stokes[0,:])

pl.show()

hsra = np.loadtxt('hsra.1d', skiprows=4)
fig, ax = pl.subplots()
ax.plot(hsra[:, 0], hsra[:, 1], label='hsra')
ax.plot(mod.atmospheres['ph1'].log_tau, mod.atmospheres['ph1'].parameters['T'], label='ph1')
ax.plot(mod.atmospheres['ph2'].log_tau, mod.atmospheres['ph2'].parameters['T'], label='ph2')

ax.legend()