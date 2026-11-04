import numpy as np
import matplotlib.pyplot as pl

if __name__ == '__main__':

    labels = ['I', 'Q', 'U', 'V']

    f, ax = pl.subplots(nrows=3, ncols=3, figsize=(15, 10), constrained_layout=True)

    mu1 = np.load('mu_1.0.npz')
    mu2 = np.load('mu_0.2.npz')

    wvl = mu1['wavelength']
    stokes1 = mu1['stokes']
    stokes2 = mu2['stokes']
    tau1 = mu1['total_tau']
    tau2 = mu2['total_tau']

    tmp = np.loadtxt('JKQ_tanausu', skiprows=1)
    j00 = tmp[:, 0]
    j20 = tmp[:, 1]
    
    PC = 2.99792458e10
    PH = 6.62606876e-27
    nbar_tanausu = (10830.0*1.e-8)**3 / (2.0*PC*PH) * j00
    omega_tanausu = j20 / j00 * np.sqrt(2.0)

    for i in range(4):
        ax.flat[i].plot(wvl, stokes1[i,:] / stokes1[0,0], label='mu=1')
        ax.flat[i].plot(wvl, stokes2[i,:] / stokes2[0,0], label='mu=0.2')
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
    ax.flat[-3].plot(J[:, 0], omega_tanausu, label='andreu')
    ax.flat[-3].set_xlabel('Optical depth')
    ax.flat[-3].set_ylabel('omega')
    ax.flat[-3].set_ylim([0.0, 0.25])
    ax.flat[-3].axhline(omega, color='r', ls='--', label='omega')
    ax.flat[-3].legend()

    ax.flat[-2].plot(wvl, tau1, label='mu=1')
    ax.flat[-2].plot(wvl, tau2, label='mu=0.2')
    ax.flat[-2].set_xlabel('Optical depth')
    ax.flat[-2].set_ylabel('Total optical depth')


    pl.savefig('test_rt.png')
    pl.show(block=True)