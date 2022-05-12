import numpy as np
import h5py

n_pixel = 100
n_lambda = 150

# Generate mock Stokes parameters full of zeros in this case (so completely useless)
stokes_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
sigma_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
los_3d = np.zeros((n_pixel,3), dtype=np.float64)
boundary_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

f = h5py.File('example.h5', 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
db_los[:] = los_3d
db_boundary[:] = boundary_3d
f.close()