import numpy as np
import hazel
import matplotlib.pyplot as pl

def read_jkq():
    f = open('JKQ_output', 'r')
    lines = f.readlines()

    j00 = []
    j20 = []

    for l in lines[1:]:
        tmp = l.split()
        tmp = [float(x) for x in tmp]
        if tmp[1] == 0.0 and tmp[2] == 0.0:
            j00.append(tmp[4])
        if tmp[1] == 2.0 and tmp[2] == 0.0:
            j20.append(tmp[4])

    return j00, j20


if __name__ == '__main__':

    # j00, j20 = read_jkq()

    # mu = 1.0
    mu = 0.2

    Bx = 10.0
    By = 10.0
    Bz = 0.0

    # Bx = 0.0
    # By = 0.0
    # Bz = 0.0

    theta = np.arccos(mu) * 180.0 / np.pi

    tmp = np.loadtxt('JKQ_tanausu', skiprows=1)
    j00 = tmp[:, 0]
    j20 = tmp[:, 1]

    PC = 2.99792458e10
    PH = 6.62606876e-27
    nbar_tanausu = (10830.0*1.e-8)**3 / (2.0*PC*PH) * j00
    omegar_tanausu = j20 / j00 * np.sqrt(2.0)

    labels = ['I', 'Q', 'U', 'V']

    mod = hazel.Model(config=None, working_mode='synthesis', verbose=3)

    mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta, 0.0, 90.0], 'Boundary condition': [1.0, 0.0, 0.0, 0.0]},0)
    mod.add_chromosphere_rt({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 8.28, 'Line': '10830', 'Wavelength': [10826, 10833], \
                        'Nslabs': 10, 'dz': 1000.0, 'nquad_photosphere': 7, 'nquad_nophotosphere': 10})
    mod.setup()

    # Vector of parameters are (Bx,By,Bz, log n, v,deltav, beta,a) and then the ff
    mod.atmospheres['ch1'].set_parameters([Bx, By, Bz, 4.34 , 0.0, 6.0, 1.0, 0.0], 1.0)

    mod.synthesize()
    

    mod2 = hazel.Model(config=None, working_mode='synthesis', verbose=3)

    mod2.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
        'LOS': [theta, 0.0, 90.0], 'Boundary condition': [1.0, 0.0, 0.0, 0.0]},0)
    mod2.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 8.28, 'Line': '10830', 'Wavelength': [10826, 10833]})
    mod2.setup()

    # Vector of parameters are (Bx,By,Bz, tau, v,deltav, beta,a) and then the ff
    mod2.atmospheres['ch1'].set_parameters([10.0, 10.0, 0.0, 2.22, 0.0, 6.0, 1.0, 0.0], 1.0)

    mod2.synthesize()

    f, ax = pl.subplots(nrows=3, ncols=3, figsize=(10, 10), constrained_layout=True)

    for i in range(4):
        ax.flat[i].plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[i,:] / mod.spectrum['spec1'].stokes[0,0], label='rt')
        ax.flat[i].plot(mod2.spectrum['spec1'].wavelength_axis, mod2.spectrum['spec1'].stokes[i,:] / mod2.spectrum['spec1'].stokes[0,0], label='hazel')
        ax.flat[i].set_xlabel('Wavelength [A]')
        ax.flat[i].set_ylabel(labels[i])


    ax.flat[0].legend()

    with open('Jbar_tensors.dat', 'r') as f:
        f.readline()
        tmp = f.readline()    
        nbar, omega = tmp.split()
        nbar = float(nbar)
        omega = float(omega)

    J = np.loadtxt('Jbar_tensors.dat', skiprows=2)
    total_tau = np.loadtxt('final_tau.dat', skiprows=1)

    ax.flat[-4].plot(J[:, 0], J[:, 2], label='hazel')
    ax.flat[-4].plot(J[:, 0], nbar_tanausu, label='andreu')
    ax.flat[-4].set_xlabel('Optical depth')
    ax.flat[-4].set_ylabel('nbar')
    ax.flat[-4].set_ylim([0.0, 0.1])
    ax.flat[-4].axhline(nbar, color='r', ls='--', label='nbar')
    ax.flat[-4].legend()

    ax.flat[-3].plot(J[:, 0], J[:, 3], label='hazel')
    ax.flat[-3].plot(J[:, 0], omegar_tanausu, label='andreu')
    ax.flat[-3].set_xlabel('Optical depth')
    ax.flat[-3].set_ylabel('omega')
    ax.flat[-3].set_ylim([0.0, 0.25])
    ax.flat[-3].axhline(omega, color='r', ls='--', label='omega')
    ax.flat[-3].legend()

    ax.flat[-2].plot(mod.spectrum['spec1'].wavelength_axis, total_tau)
    ax.flat[-2].set_xlabel('Optical depth')
    ax.flat[-2].set_ylabel('Total optical depth')

    np.savez(f'mu_{mu}_Bx_{Bx}_By_{By}_Bz_{Bz}.npz', 
             wavelength=mod.spectrum['spec1'].wavelength_axis, 
             stokes=mod.spectrum['spec1'].stokes / mod.spectrum['spec1'].stokes[0,0], 
             total_tau=total_tau,
             tau=J[:, 0],
             nbar=J[:, 2],
             omega=J[:, 3])


    # pl.savefig('test_rt.png')
    pl.show(block=True)

    # REVISAR LINEA 137 rt_coef.f90
