import hazel
import matplotlib.pyplot as pl


mod = hazel.Model(config=None, working_mode='synthesis', verbose=3)

mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
    'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]},0)
mod.add_chromosphere_rt({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 8.28, 'Line': '10830', 'Wavelength': [10826, 10833], \
                      'Nslabs': 10, 'dz': 10000.0, 'nquad_photosphere': 3, 'nquad_nophotosphere': 7})
mod.setup()

# Vector of parameters are (Bx,By,Bz, log n, v,deltav, beta,a) and then the ff
mod.atmospheres['ch1'].set_parameters([10.0, 10.0, 0.0, 3.7 , 0.0, 8.0, 1.0, 0.0], 1.0)

mod.synthesize()


mod2 = hazel.Model(config=None, working_mode='synthesis', verbose=3)

mod2.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
    'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]},0)
mod2.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 8.28, 'Line': '10830', 'Wavelength': [10826, 10833]})
mod2.setup()

# Vector of parameters are (Bx,By,Bz, tau, v,deltav, beta,a) and then the ff
mod2.atmospheres['ch1'].set_parameters([10.0, 10.0, 0.0, 3.5 , 0.0, 8.0, 1.0, 0.0], 1.0)

mod2.synthesize()

f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()

for i in range(4):
    ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
    ax[i].plot(mod2.spectrum['spec1'].stokes[i,:])

pl.show(block=True)