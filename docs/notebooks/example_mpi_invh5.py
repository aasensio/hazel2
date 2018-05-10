import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop


iterator = hazel.iterator(use_mpi=True)
rank = iterator.get_rank()

if (rank == 0):    
    mod = hazel.Model('conf_mpi_invh5.ini', working_mode='inversion')
    iterator.use_model(model=mod)
else:
    iterator.use_model()

iterator.run_all_pixels()