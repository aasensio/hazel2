import numpy as np
import matplotlib.pyplot as pl
import hazel

mod = hazel.Model('conf_spot.ini', working_mode='inversion', verbose=3)
mod.read_observation()
mod.open_output()
mod.invert()
mod.write_output()
mod.close_output()