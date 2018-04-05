import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
# from ipdb import set_trace as stop


# Test a single inversion in non-iterator mode
mod = hazel.Model('conf_single.ini', working_mode='synthesis')
mod.synthesize()

assert mod.spectrum['spec1'].stokes.shape == (4,150), "incorrect dimensions in synthesis"
# f, ax = pl.subplots(nrows=2, ncols=2)
# ax = ax.flatten()
# for i in range(4):
#     ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
# pl.show()
# stop()