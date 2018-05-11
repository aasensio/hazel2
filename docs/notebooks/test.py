import numpy as np
import matplotlib.pyplot as pl
import hazel

mod = hazel.Model('conf_single.ini', working_mode='synthesis')
mod.open_output()
mod.synthesize()
mod.write_output()
mod.close_output()
