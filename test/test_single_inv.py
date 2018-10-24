import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py

def test_single_inv():
# Test a single inversion in non-iterator mode without randomization
    mod = hazel.Model('test/configurations/conf_single.ini', working_mode='inversion', verbose=2)
    mod.read_observation()
    mod.open_output()

    mod.invert()
    mod.write_output()

    mod.close_output()

    f = h5py.File('output.h5')

    assert f['spec1']['stokes'].shape == (1,1,2,4,150), "incorrect dimensions in inversion"

    f.close()