import numpy as np
import hazel

def test_single_syn():
    # Test a single inversion in non-iterator mode
    mod = hazel.Model('conf_single.ini', working_mode='synthesis')
    mod.synthesize()

    assert mod.spectrum['spec1'].stokes.shape == (4,150), "incorrect dimensions in synthesis"