import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py

# Test iterator with a single observation in synthesis
iterator = hazel.Iterator(use_mpi=False)
rank = iterator.get_rank()
mod = hazel.Model('configurations/conf_nonmpi_inv1d.ini', working_mode='inversion', verbose=2, randomization=2)
iterator.use_model(model=mod)
iterator.run_all_pixels()

if (rank == 0):
    fig, ax = pl.subplots(nrows=2, ncols=2)
    ax = ax.flatten()

    f = h5py.File('output.h5', 'r')
    obs = np.loadtxt('observations/10830_stokes.1d', skiprows=7)

    for i in range(4):
        ax[i].plot(obs[:,i])
        ax[i].plot(f['spec1']['stokes'][0,0,0,i,:])
        ax[i].plot(f['spec1']['stokes'][0,1,0,i,:])

    pl.show()
    pl.pause(0.001)

    f.close()

    input("Press [enter] to continue.")
