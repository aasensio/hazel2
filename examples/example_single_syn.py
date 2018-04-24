import numpy as np
import matplotlib.pyplot as pl
import hazel

# Test a single inversion in non-iterator mode
mod = hazel.Model('conf_single.ini', working_mode='synthesis')
mod.open_output()
mod.synthesize()
mod.write_output()
mod.close_output()

f, ax = pl.subplots()
ax.plot(mod.spectrum['spec1'].stokes[0,:])
pl.show()
