import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py

tmp = hazel.tools.File_photosphere(mode='single')
tmp.set_default(n_pixel=1)
tmp.save('photospheres/init')

tmp = hazel.tools.File_chromosphere(mode='single')
tmp.set_default(n_pixel=1)
tmp.save('chromospheres/init')

mod = hazel.Model('conf_spot.ini', working_mode='inversion', verbose=3)
mod.read_observation()
mod.open_output()
mod.invert()
mod.write_output()
mod.close_output()

f = h5py.File('output.h5', 'r')

stokes = np.loadtxt('10830_spot.1d', skiprows=7).T

print('(npix,nrand,ncycle,nstokes,nlambda) -> {0}'.format(f['spec1']['stokes'].shape))

fig, ax = pl.subplots(nrows=2, ncols=2, figsize=(10,10))
ax = ax.flatten()
for i in range(4):
    ax[i].plot(f['spec1']['wavelength'][:] - 10830, stokes[i,:])
    for j in range(2):
        ax[i].plot(f['spec1']['wavelength'][:] - 10830, f['spec1']['stokes'][0,0,j,i,:])

for i in range(4):
    ax[i].set_xlabel('Wavelength - 10830[$\AA$]')    
    ax[i].set_xlim([-7,3])

pl.tight_layout()

f.close()
