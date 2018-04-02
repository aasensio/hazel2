import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

# Test iterator with a single observation in synthesis
iterator = hazel.iterator(use_mpi=False)
rank = iterator.get_rank()
mod = hazel.Model('conf_nonmpi_syn1d.ini')
iterator.use_model(model=mod)
iterator.run_all_pixels()

if (rank == 0):
    f = h5py.File('output.h5', 'r')
    pl.plot(f['spec1'][:,0,:].T/f['spec1'][:,0,-2][:,None].T)
    pl.show()
    pl.pause(0.001)

    input("Press [enter] to continue.")