import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

# Test iterator with a single observation in synthesis
iterator = hazel.iterator(use_mpi=False)
rank = iterator.get_rank()
mod = hazel.Model('conf_nonmpi_inv1d.ini', working_mode='inversion', verbose=3)
iterator.use_model(model=mod)
iterator.run_all_pixels()

if (rank == 0):
    fig, ax = pl.subplots(nrows=2, ncols=2)
    ax = ax.flatten()

    f = h5py.File('output.h5', 'r')

    for i in range(4):
        ax[i].plot(f['spec1'][0,0,i,:])

    pl.show()
    pl.pause(0.001)

    input("Press [enter] to continue.")
