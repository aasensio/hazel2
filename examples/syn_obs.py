import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop


def read_photosphere(file):
    f = open(file, 'r')
    f.readline()
    ff = float(f.readline())
    f.close()
    return np.loadtxt(file, skiprows=4), ff

def read_1d(file):
    f = open(file, 'r')
    f.readline()
    los = np.array(f.readline().split()).astype('float64')
    f.readline()
    f.readline()
    boundary = np.array(f.readline().split()).astype('float64')
    f.readline()
    f.readline()
    tmp = f.readlines()
    f.close()
    n_lambda = len(tmp)
    stokes = np.zeros((4,n_lambda))
    noise = np.zeros((4,n_lambda))
    for i, l in enumerate(tmp):
        t = np.array(l.split()).astype('float64')
        stokes[:,i] = t[0:4]
        noise[:,i] = t[4:]

    mu = np.cos(los[0] * np.pi / 180.0)

    return stokes, noise, los, mu, boundary

if __name__ == '__main__':

    #-------------------------------------------
    # Synthesize a sample profile with some default configuration given in conf_syn.ini
    # The observed files will be used later as observations for the inversions
    #-------------------------------------------
    mod = hazel.Model('configurations/conf_syn.ini')
    mod.synthesize()


    noise = 1e-4 * np.ones((150,4))
    spec = mod.spectrum['spec1'].stokes.T + 1e-4 * np.random.randn(150,4)

    # Generate wavelength axis
    np.savetxt('observations/10830.wavelength', mod.spectrum['spec1'].wavelength_axis, header='lambda')

    # Generate output Stokes parameters
    f = open('observations/10830_stokes.1d', 'w')
    f.write('# LOS theta_LOS, phi_LOS, gamma_LOS\n')
    f.write('0.0 0.0 90.0\n')
    f.write('\n')
    f.write('# Boundary condition I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)\n')
    f.write('1.0 0.0 0.0 0.0\n')
    f.write('\n')
    f.write('# SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV\n')
    for i in range(150):
        f.write('{0} {1}\n'.format('  '.join(spec[i,:].astype('str')), '  '.join(noise[i,:].astype('str'))))
    f.close()

    # Files with wavelength weights
    f = open('observations/10830.weights', 'w')
    f.write('# WeightI WeightQ WeightU WeightV\n')
    for i in range(150):
        f.write('1.0    1.0   1.0   1.0\n')
    f.close()


    #-------------------------------------------
    # Now we generate HDF5 files with the same pixel repeated several times
    #-------------------------------------------
    n_pixel = 10

    model, ff = read_photosphere('photospheres/model_photosphere.1d')

    nz, ncol = model.shape
    model_3d = np.zeros((n_pixel,nz,ncol), dtype=np.float64)
    ff_3d = np.zeros((n_pixel,), dtype=np.float64)

    for i in range(n_pixel):
        model[:,1] += 100.0
        model_3d[i,:,:] = model
        ff_3d[i] = ff
        
    f = h5py.File('photospheres/model_photosphere.h5', 'w')
    db_model = f.create_dataset('model', model_3d.shape, dtype=np.float64)
    db_ff = f.create_dataset('ff', ff_3d.shape, dtype=np.float64)
    db_model[:] = model_3d
    db_ff[:] = ff_3d
    f.close()

    # Observations    
    stokes, sigma, los, mu, boundary = read_1d('observations/10830_stokes.1d')

    n_lambda, _ = stokes.T.shape

    stokes_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
    sigma_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
    los_3d = np.zeros((n_pixel,3), dtype=np.float64)
    boundary_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

    for i in range(n_pixel):
        stokes_3d[i,:,:] = stokes.T
        sigma_3d[i,:,:] = sigma.T
        los_3d[i,:] = np.array([0.0,0.0,90.0])
        boundary_3d[i,:,:] = np.repeat(np.atleast_2d(boundary), n_lambda, axis=0)

    f = h5py.File('observations/10830_stokes.h5', 'w')
    db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
    db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
    db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
    db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
    db_stokes[:] = stokes_3d
    db_sigma[:] = sigma_3d
    db_los[:] = los_3d
    db_boundary[:] = boundary_3d
    f.close()


    model = np.loadtxt('chromospheres/model_chromosphere.1d', skiprows=1)

    model_3d = np.zeros((n_pixel,8), dtype=np.float64)
    ff_3d = np.zeros((n_pixel,), dtype=np.float64)

    for i in range(n_pixel):    
        model_3d[i,:] = model[0:-1]
        ff_3d[i] = model[-1]

    f = h5py.File('photospheres/model_chromosphere.h5', 'w')
    db_model = f.create_dataset('model', model_3d.shape, dtype=np.float64)
    db_ff = f.create_dataset('ff', ff_3d.shape, dtype=np.float64)
    db_model[:] = model_3d
    db_ff[:] = ff_3d
    f.close()


