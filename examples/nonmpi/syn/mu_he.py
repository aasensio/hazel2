import hazel
import numpy as np
import matplotlib.pyplot as pl

mod = hazel.Model(working_mode='synthesis', verbose=3)
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1', 'LOS': [0.0, 0.0, 90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]}, index=0)
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10833], 'Reference frame': 'vertical'})
mod.setup()

theta = [0.0, 30.0, 60.0]

mod.atmospheres['ch1'].set_parameters([0.0, 0.0, 0.0,1.0,0.0,8.0,1.0,0.0],1.0)

fig, ax = pl.subplots()

for i in range(len(theta)):
    mod.spectrum['spec1'].set_los(np.array([theta[i], 0.0, 90.0]))
    mod.synthesize()
    ax.plot(mod.spectrum['spec1'].stokes[0, :] / mod.spectrum['spec1'].stokes[0, 0], label=str(theta[i])+' deg')

ax.legend()
