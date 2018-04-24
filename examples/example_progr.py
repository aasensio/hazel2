import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
#from ipdb import set_trace as stop


# Test a single synthesis in programmatic mode
mod = hazel.Model(working_mode='synthesis')
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1', 'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10833]})
mod.setup()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()

for j in range(5):
    mod.atmospheres['ch1'].set_parameters([0.0,0.0,100.0*j,1.0,0.0,8.0,1.0,0.0],1.0)
    mod.synthesize()

    for i in range(4):
        ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
pl.show()
