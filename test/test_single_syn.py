import numpy as np
import matplotlib.pyplot as pl
import hazel

def test_single_syn():
# Test a single inversion in non-iterator mode
    mod = hazel.Model('test/configurations/conf_single.ini', working_mode='synthesis', verbose=3)
    mod.open_output()
    mod.synthesize()
    mod.write_output()
    mod.close_output()

    assert mod.spectrum['spec1'].stokes.shape == (4,150), "incorrect dimensions in synthesis"    