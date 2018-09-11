import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop


iterator = hazel.Iterator(use_mpi=True)
rank = iterator.get_rank()
mod = hazel.Model('conf_mpi_synh5.ini', working_mode='synthesis', rank=rank)
iterator.use_model(model=mod)
iterator.run_all_pixels()