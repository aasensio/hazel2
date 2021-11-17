import numpy as np
import matplotlib.pyplot as pl
import hazel

tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1, default='falc82')
tmp.save('../../photospheres/falc')


# Test a single inversion in non-iterator mode
mod = hazel.Model('../../configurations/conf_caii.ini', working_mode='synthesis', verbose=3, root='../../')
mod.synthesize(nlte=False)

stokes_lte = np.copy(mod.spectrum['spec1'].stokes[0,:])

mod.synthesize(nlte=True)
stokes_nlte = mod.spectrum['spec1'].stokes[0,:]

fig, ax = pl.subplots(nrows=1, ncols=2, figsize=(12, 6))
ax[0].plot(mod.atmospheres['ph1'].log_tau, np.log10(mod.atmospheres['ph1'].departure[0, 0, :]))
ax[0].plot(mod.atmospheres['ph1'].log_tau, np.log10(mod.atmospheres['ph1'].departure[1, 0, :]))

ax[0].set_ylim([-1,9])

ax[1].plot(mod.spectrum['spec1'].wavelength_axis, stokes_lte, label='LTE')
ax[1].plot(mod.spectrum['spec1'].wavelength_axis, stokes_nlte, label='NLTE')

ax[1].legend()

pl.show()