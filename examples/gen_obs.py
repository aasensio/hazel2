import numpy as np
import hazel
import h5py

mod = hazel.Model('configurations/conf_syn.ini', working_mode='synthesis')
mod.open_output()
mod.synthesize()
mod.write_output()
mod.close_output()


f = h5py.File('output.h5', 'r')

stokes = f['spec1']['stokes']

tmp = hazel.tools.File_observation(mode='single')
tmp.set_size(n_lambda=150, n_pixel=1)

tmp.obs['stokes'][0, :, :] = stokes[0, 0, :, :].T
tmp.obs['sigma'][0, :, :] = 1e-3
tmp.obs['los'][0, :] = np.array([0.0, 0.0, 90.0])
tmp.obs['boundary'][0, :] = np.array([1.0, 0.0, 0.0, 0.0])


tmp.save('test')


f.close()