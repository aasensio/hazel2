import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

version = 'serial_invert'


# Single invert
if (version == 'single_invert'):
    mod = hazel.Model('conf_single_invert.ini')
    mod.invert()

    # pl.plot(mod.spectrum['spec1'].stokes[0,:])
    # pl.show()
    # pl.pause(0.001)

    # np.savetxt('10830_stokes.1d', mod.spectrum['spec1'].stokes.T + 1e-4 * np.random.randn(150,4))
    # input("Press [enter] to continue.")
    

# Serial invert
if (version == 'serial_invert'):
    iterator = hazel.iterator(use_mpi=False)
    rank = iterator.get_rank()
    mod = hazel.Model('conf_single_invert.ini')
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

# Serial
if (version == 'serial'):
    iterator = hazel.iterator(use_mpi=False)
    rank = iterator.get_rank()
    mod = hazel.Model('conf.ini')
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

    if (rank == 0):
        f = h5py.File('output.h5', 'r')
        pl.plot(f['spec1'][:,0,:].T/f['spec1'][:,0,-2][:,None].T)
        pl.show()
        pl.pause(0.001)

        input("Press [enter] to continue.")

# Parallel
if (version == 'parallel'):
    iterator = hazel.iterator(use_mpi=True)
    rank = iterator.get_rank()

    if (rank == 0):    
        mod = hazel.Model('conf.ini')
        iterator.use_model(model=mod)
    else:
        iterator.use_model()

    iterator.run_all_pixels()

    # if (rank == 0):
    #     f = h5py.File('output.h5', 'r')
    #     pl.plot(f['spec1'][:,0,:].T/f['spec1'][:,0,-2][:,None].T)
    #     pl.show()
    #     pl.pause(0.001)

    #     input("Press [enter] to continue.")