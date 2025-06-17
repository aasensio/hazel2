import hazel
import numpy as np
import matplotlib.pyplot as pl

mod = hazel.Model(working_mode='synthesis', verbose=3)
mod.add_spectral({'Name': 'spec1', 'Wavelength': [6300.8921, 6303.2671, 112], 'topology': 'ph1',
    'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]}, index=0)
mod.add_photosphere({'Name': 'ph1', 'Spectral region': 'spec1', 'Spectral lines': [200, 201],
                        'Wavelength': [6300.8921, 6303.2671]})
mod.setup()

model = np.loadtxt('../../photospheres/model_photosphere.1d', skiprows=4)

model[:, 2] = -1.0

mod.atmospheres['ph1'].set_parameters(model, 1.0, 0.0)

theta = [0.0, 30.0, 60.0]

fig, ax = pl.subplots()

for i in range(len(theta)):
    mod.spectrum['spec1'].set_los([theta[i], 0.0, 90.0])
    mod.synthesize()
    ax.plot(mod.spectrum['spec1'].stokes[0, :] / mod.spectrum['spec1'].stokes[0, 0])
    # ax.plot(mod.spectrum['spec1'].stokes[0, :])
