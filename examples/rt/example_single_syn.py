import numpy as np
import matplotlib.pyplot as pl
import hazel

# Test a single inversion in non-iterator mode
mod_rt = hazel.Model('configurations/conf_single_rt.ini', working_mode='synthesis', verbose=3, root='./')
mod_rt.open_output()
mod_rt.synthesize()
mod_rt.write_output()
mod_rt.close_output()

mod = hazel.Model('configurations/conf_single.ini', working_mode='synthesis', verbose=3, root='./')
mod.open_output()
mod.synthesize()
mod.write_output()
mod.close_output()

fig, ax = pl.subplots(nrows=2, ncols=2, figsize=(10,10))
for i in range(4):
    ax.flat[i].plot(mod.spectrum['spec1'].stokes[i,:])
    ax.flat[i].plot(mod_rt.spectrum['spec1'].stokes[i,:])
pl.show()
