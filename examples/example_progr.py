import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
#from ipdb import set_trace as stop

reference = 'line-of-sight'
# reference = 'vertical'

theta = 45.0

if (reference == 'line-of-sight'):    
    thetaB = 0.0
    phiB = 0.0

if (reference == 'vertical'):
    thetaB = 45.0
    phiB = 0.0    

# Test a single synthesis in programmatic mode
mod = hazel.Model(working_mode='synthesis')
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1', 'LOS': [theta,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10833], 'Reference frame': reference})
mod.setup()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()

for j in range(5):

    B = 100.0 * j + 10    

    Bx = B * np.sin(thetaB*np.pi/180) * np.cos(phiB*np.pi/180.)
    By = B * np.sin(thetaB*np.pi/180) * np.sin(phiB*np.pi/180.)
    Bz = B * np.cos(thetaB*np.pi/180)
    
    mod.atmospheres['ch1'].set_parameters([Bx, By, Bz,1.0,0.0,8.0,1.0,0.0],1.0)
    mod.synthesize()

    for i in range(4):
        ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
pl.show()
