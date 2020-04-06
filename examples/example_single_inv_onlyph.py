import numpy as np
import matplotlib.pyplot as pl
import hazel

# Test a single inversion in non-iterator mode
mod = hazel.Model('configurations/conf_single_onlyph.ini', working_mode='inversion', verbose=4)
mod.read_observation()
mod.open_output()

mod.invert()
mod.write_output()

mod.close_output()

f, ax = pl.subplots()
ax.plot(mod.spectrum['spec1'].obs[0,:])
ax.plot(mod.spectrum['spec1'].stokes[0,:])
pl.show()
