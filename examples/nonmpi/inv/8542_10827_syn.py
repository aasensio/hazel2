import numpy as np
import matplotlib.pyplot as pl
import hazel

tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1, default='falc82')
tmp.save('falc')

# Test a single inversion in non-iterator mode
mod = hazel.Model('conf_8542_10827_syn.ini', working_mode='synthesis', verbose=4, root='./')
mod.synthesize()

fig, ax = pl.subplots(nrows=1, ncols=2, figsize=(12, 6))

ax[0].plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[0,:])
ax[1].plot(mod.spectrum['spec2'].wavelength_axis, mod.spectrum['spec2'].stokes[0,:])


pl.show()

tmp = hazel.tools.File_observation(mode='single')
tmp.set_size(n_lambda=150, n_pixel=1)

tmp.obs['stokes'][0, :, :] = mod.spectrum['spec1'].stokes.T
tmp.obs['sigma'][0, :, :] = 1e-3
tmp.obs['wavelength'] = mod.spectrum['spec1'].wavelength_axis
tmp.obs['boundary'][0, :, 0] = 1.0
tmp.obs['los'][0, :] = np.array([0.0, 0.0, 90.0])
tmp.obs['weights'] = np.ones((150, 4))
tmp.obs['mask'] = np.ones((1,))
tmp.save('8542')

tmp.obs['stokes'][0, :, :] = mod.spectrum['spec2'].stokes.T
tmp.obs['sigma'][0, :, :] = 1e-3
tmp.obs['wavelength'] = mod.spectrum['spec2'].wavelength_axis
tmp.obs['boundary'][0, :, 0] = 1.0
tmp.obs['los'][0, :] = np.array([0.0, 0.0, 90.0])
tmp.obs['weights'] = np.ones((150, 4))
tmp.obs['mask'] = np.ones((1,))
tmp.save('10827')