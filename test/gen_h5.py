import numpy as np
import h5py

def read_photosphere(file):
    f = open(file, 'r')
    f.readline()
    ff = float(f.readline())
    f.close()
    return np.loadtxt(file, skiprows=4), ff

n_pixel = 50

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
n_pixel = 2
obs = np.loadtxt('observations/10830_stokes.1d', skiprows=1)
stokes = obs[:,0:4]
sigma = obs[:,4:]

n_lambda, _ = stokes.shape

stokes_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
sigma_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

for i in range(n_pixel):
    stokes_3d[i,:,:] = stokes
    sigma_3d[i,:,:] = sigma

f = h5py.File('observations/10830_stokes.h5', 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
f.close()


# model = np.loadtxt('chromospheres/model_chromosphere.1d', skiprows=1)

# model_3d = np.zeros((n_pixel,8), dtype=np.float64)
# for i in range(n_pixel):
#     model_3d[i,:] = model

# f = h5py.File('chromospheres/model_chromosphere.h5', 'w')
# db = f.create_dataset('model', model_3d.shape, dtype=np.float64)
# db[:] = model_3d
# f.close()

# model = np.loadtxt('chromospheres/model_chromosphere2.1d', skiprows=1)

# model_3d = np.zeros((n_pixel,8), dtype=np.float64)
# for i in range(n_pixel):
#     model[5] += 2.0
#     model_3d[i,:] = model
    
# f = h5py.File('chromospheres/model_chromosphere2.h5', 'w')
# db = f.create_dataset('model', model_3d.shape, dtype=np.float64)
# db[:] = model_3d
# f.close()