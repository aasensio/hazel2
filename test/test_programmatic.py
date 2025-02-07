import numpy as np
import hazel


def test_programmatic():
    # Test a single inversion in non-iterator mode
    mod = hazel.Model(working_mode='synthesis')
    mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1', 'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]}i, 0)
    mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10833]})
    mod.setup()
    mod.atmospheres['ch1'].set_parameters([0.0,0.0,100.0,1.0,0.0,8.0,1.0,0.0],1.0)
    mod.synthesize()

    assert mod.spectrum['spec1'].stokes.shape == (4,150), "incorrect dimensions in synthesis"
