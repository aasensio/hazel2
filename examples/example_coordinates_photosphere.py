import hazel
import numpy as np
import matplotlib.pyplot as pl

# Define LOS theta angle
theta_los = 25.0

# Atmosphere with vertical/cartesian
mod = hazel.Model(working_mode='synthesis', verbose=3)
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10824, 10830, 150], 'topology': 'ph1',
    'LOS': [theta_los,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
mod.add_photosphere({'Name': 'ph1', 'Spectral region': 'spec1', 'Spectral lines': [300,],
    'Wavelength': [10824, 10833], 'Reference atmospheric model': 'photospheres/model_photosphere.1d'})
mod.setup()

mod.synthesize()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()

for i in range(4):
    ax[i].plot(mod.spectrum['spec1'].stokes[i,:])


model = mod.atmospheres['ph1'].get_parameters()

model[:,5] = 0.0
model[:,6] = 500.0
mod.atmospheres['ph1'].set_parameters(model, 1.0, 0.0)

mod.synthesize()

for i in range(4):
    ax[i].plot(mod.spectrum['spec1'].stokes[i,:])

pl.show()