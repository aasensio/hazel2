import numpy as np
import matplotlib.pyplot as pl
import hazel

tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1, default='hsra')
tmp.save('../../photospheres/falc')


# Test a single inversion in non-iterator mode
mod = hazel.Model('../../configurations/conf_caii.ini', working_mode='synthesis', verbose=3, root='../../')
mod.open_output()
mod.synthesize(nlte=False)
mod.write_output()
mod.close_output()

f, ax = pl.subplots()
ax.plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[0,:])


mod.synthesize(nlte=True)
ax.plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[0,:])

pl.show()
