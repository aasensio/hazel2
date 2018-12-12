import numpy as np
import matplotlib.pyplot as pl
import hazel
#from ipdb import set_trace as stop

reference = 'line-of-sight'
# reference = 'vertical'

theta = 85.0

if (reference == 'line-of-sight'):    
    thetaB = 0.0
    phiB = 0.0

if (reference == 'vertical'):
    thetaB = 45.0
    phiB = 0.0    

pl.close('all')

# Test a photosphere
mod = hazel.Model(working_mode='synthesis')
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10835, 150], 'topology': 'ph1', 'LOS': [theta,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
mod.add_photosphere({'Name': 'ph1', 'Spectral region': 'spec1', 'Spectral lines': [300], 'Wavelength': [10826, 10835], 'Reference atmospheric model': 'photospheres/model_photosphere.1d'})
mod.setup()
mod.synthesize()

f, ax = pl.subplots()
ax.plot(mod.spectrum['spec1'].stokes[0,:])


mod = hazel.Model(working_mode='synthesis')
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10835, 150], 'topology': 'ch1', 'LOS': [theta,0.0,90.0], 'Boundary condition': [0.5,0.0,0.0,0.0]})
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10835], 'Reference frame': reference})
mod.setup()


mod.atmospheres['ch1'].set_parameters([0.0, 0.0, 0.0,1.0,0.0,8.0,1.0,0.0],1.0)
mod.synthesize()

ax.plot(mod.spectrum['spec1'].stokes[0,:])

mod = hazel.Model(working_mode='synthesis')
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10835, 150], 'topology': 'ph1->ch1', 'LOS': [theta,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
mod.add_photosphere({'Name': 'ph1', 'Spectral region': 'spec1', 'Spectral lines': [300], 'Wavelength': [10826, 10835], 'Reference atmospheric model': 'photospheres/model_photosphere.1d'})
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10835], 'Reference frame': reference})
mod.setup()

mod.atmospheres['ch1'].set_parameters([0.0, 0.0, 0.0,1.0,0.0,8.0,1.0,0.0],1.0)
mod.synthesize()


ax.plot(mod.spectrum['spec1'].stokes[0,:])
pl.show()