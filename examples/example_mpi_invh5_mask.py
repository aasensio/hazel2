import numpy as np
import matplotlib.pyplot as pl
import hazel

iterator = hazel.Iterator(use_mpi=True)
rank = iterator.get_rank()
mod = hazel.Model('configurations/conf_mpi_invh5_mask.ini', working_mode='inversion', rank=rank, verbose=3)
iterator.use_model(model=mod)
iterator.run_all_pixels()
