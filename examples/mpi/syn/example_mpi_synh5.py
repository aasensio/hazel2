import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py


iterator = hazel.Iterator(use_mpi=True)
rank = iterator.get_rank()
mod = hazel.Model('../../configurations/conf_mpi_synh5.ini', working_mode='synthesis', rank=rank, root='../../')
iterator.use_model(model=mod)
iterator.run_all_pixels()
