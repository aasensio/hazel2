import numpy as np
import h5py
from ipdb import set_trace as stop

__all__ = ['File_observation', 'File_photosphere', 'File_chromosphere']

class File_observation(object):
    """
    Class that defines an observation. This can be used to easily save observations in the appropriate format
    """
    def __init__(self, file=None, mode='single'):

        if (file is not None):
            self.obs = self.read(file)

        self.mode = mode
        
        self.obs = {'stokes': None, 'sigma': None, 'los': None, 'boundary': None, 'wavelength': None, 'weights': None}

    def set_size(self, n_lambda, n_pixel=1):
        """
        Set the number of wavelengths and number of pixels of the current observation

        Parameters
        ----------
        n_lambda : int
            Number of wavelength points

        n_pixel : int (optional, equal to 1 as default)
            Number of pixels of the output

        Returns
        -------
        None
        """     

        if (self.mode == 'single' and n_pixel > 1):
            raise Exception("Single pixel models cannot contain more than one pixel")

        self.n_pixel = n_pixel
        self.n_lambda = n_lambda

        self.obs['stokes'] = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
        self.obs['sigma'] = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
        self.obs['los'] = np.zeros((n_pixel,3), dtype=np.float64)
        self.obs['boundary'] = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

        self.obs['mask'] = np.zeros((n_pixel,), dtype=np.int8)

        self.obs['weights'] = np.zeros((n_lambda,4), dtype=np.float64)
        self.obs['wavelength'] = np.zeros((n_lambda), dtype=np.float64)
        
    def save(self, file):
        """
        Save the curent observation

        Parameters
        ----------
        file : str
            Name of the output files. Extensions will be added to it
        
        Returns
        -------
        None
        """
        
        # Save wavelength
        print("Saving wavelength file : {0}.wavelength".format(file))
        np.savetxt('{0}.wavelength'.format(file), self.obs['wavelength'], header='lambda')

        # Save weights
        print("Saving weights file : {0}.weights".format(file))
        f = open('{0}.weights'.format(file), 'w')
        f.write('# WeightI WeightQ WeightU WeightV\n')
        for i in range(self.n_lambda):
            f.write('{0}   {1}    {2}    {3}\n'.format(self.obs['weights'][i,0], self.obs['weights'][i,1], self.obs['weights'][i,2], self.obs['weights'][i,3]))
        f.close()

        if (self.mode == 'single'):
            print("Saving 1D Stokes file : {0}.1d".format(file))
            f = open('{0}.1d'.format(file), 'w')
            f.write('# LOS theta_LOS, phi_LOS, gamma_LOS\n')
            f.write('{0}   {1}   {2}\n'.format(self.obs['los'][0,0], self.obs['los'][0,1], self.obs['los'][0,2]))
            f.write('\n')
            f.write('# Boundary condition I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)\n')
            f.write('{0}   {1}   {2}    {3}\n'.format(self.obs['boundary'][0,0,0], self.obs['boundary'][0,0,1], self.obs['boundary'][0,0,2], self.obs['boundary'][0,0,3]))
            f.write('\n')
            f.write('# SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV\n')
            tmp = np.hstack([np.squeeze(self.obs['stokes']), np.squeeze(self.obs['sigma'])])
            np.savetxt(f, tmp)
            f.close()

        if (self.mode == 'multi'):
            print("Saving 3D Stokes file : {0}.h5".format(file))
            f = h5py.File('{0}.h5'.format(file), 'w')
            db_stokes = f.create_dataset('stokes', self.obs['stokes'].shape, dtype=np.float64)
            db_sigma = f.create_dataset('sigma', self.obs['sigma'].shape, dtype=np.float64)
            db_los = f.create_dataset('LOS', self.obs['los'].shape, dtype=np.float64)
            db_boundary = f.create_dataset('boundary', self.obs['boundary'].shape, dtype=np.float64)
            db_stokes[:] = self.obs['stokes']
            db_sigma[:] = self.obs['sigma']
            db_los[:] = self.obs['los']
            db_boundary[:] = self.obs['boundary']
            f.close()

            print("Saving 3D mask file : {0}.mask".format(file))
            f = h5py.File('{0}.mask'.format(file), 'w')
            db_mask = f.create_dataset('mask', self.obs['mask'].shape, dtype=np.int8)            
            db_mask[:] = self.obs['mask']
            f.close()


class File_photosphere(object):
    """
    Class that defines a model photosphere and can be used to easily save observations
    """
    def __init__(self, file=None, mode='single'):

        if (file is not None):
            self.model = self.read(file)

        self.mode = mode
        
        self.model = {'model': None, 'ff': None}

    def set_size(self, nz, n_pixel=1):
        """
        Set the number of depth points and number of pixels of the current atmosphere

        Parameters
        ----------
        nz : int
            Number of depth points of the atmosphere

        n_pixel : int (optional, equal to 1 as default)
            Number of pixels of the output

        Returns
        -------
        None
        """

        if (self.mode == 'single' and n_pixel > 1):
            raise Exception("Single pixel models cannot contain more than one pixel")

        self.n_pixel = n_pixel
        self.n_lambda = n_lambda

        self.model['model'] = np.zeros((n_pixel,nz,8), dtype=np.float64)
        self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)

    def set_default(self, n_pixel=1, default='hsra'):
        """
        Set the atmosphere to one of the default ones available in the code

        Parameters
        ----------        
        n_pixel : int (optional, equal to 1 as default)
            Number of pixels of the output

        default : str ('hsra' -> Harvard-Smithsonian Reference Atmosphere)

        Returns
        -------
        None
        """
        if (self.mode == 'single' and n_pixel > 1):
            raise Exception("Single pixel models cannot contain more than one pixel")

        if (default == 'hsra'):
            print("Setting HSRA photosphere")
            path = str(__file__).split('/')
            filename = '/'.join(path[0:-1])+'/data/hsra.mod'
            f = open(filename, 'r')
            f.readline()
            ff = float(f.readline())
            f.close()
            model = np.loadtxt(filename, skiprows=4)

            nz = model.shape[0]

            self.model['model'] = np.zeros((n_pixel,nz,8), dtype=np.float64)
            self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)

            self.model['model'][:] = model[None,:,:]
            self.model['ff'][:] = ff
        
    def save(self, file, default=None):
        """
        Save the curent model

        Parameters
        ----------
        file : str
            Name of the output files. Extensions will be added to it        

        Returns
        -------
        None
        """
        
        if (self.mode == 'single'):
            print("Saving photospheric 1D model : {0}.1d".format(file))
            f = open('{0}.1d'.format(file), 'w')
            f.write('ff\n')
            f.write('{0}\n'.format(self.model['ff'][0]))
            f.write('\n')
            f.write('  logtau     T        Pe           vmic        v            Bx           By         Bz\n')            
            
            np.savetxt(f, np.squeeze(self.model['model']))
            f.close()

        if (self.mode == 'multi'):
            print("Saving photospheric 3D model : {0}.h5".format(file))
            f = h5py.File('{0}.h5'.format(file), 'w')
            db_model = f.create_dataset('model', self.model['model'].shape, dtype=np.float64)
            db_ff = f.create_dataset('ff', self.model['ff'].shape, dtype=np.float64)
            db_model[:] = self.model['model']
            db_ff[:] = self.model['ff']
            f.close()


class File_chromosphere(object):
    """
    Class that defines a model atmosphere and can be used to easily save observations
    """
    def __init__(self, file=None, mode='single'):

        if (file is not None):
            self.model = self.read(file)

        self.mode = mode
        
        self.model = {'model': None, 'ff': None}

    def set_size(self, n_pixel=1):
        """
        Set the number of pixels of the current atmosphere

        Parameters
        ----------        
        n_pixel : int (optional, equal to 1 as default)
            Number of pixels of the output

        Returns
        -------
        None
        """

        if (self.mode == 'single' and n_pixel > 1):
            raise Exception("Single pixel models cannot contain more than one pixel")

        self.n_pixel = n_pixel
        self.n_lambda = n_lambda

        self.model['model'] = np.zeros((n_pixel,8), dtype=np.float64)
        self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)

    def set_default(self, n_pixel=1, default='disk'):
        """
        Set the atmosphere to one of the default ones available in the code

        Parameters
        ----------        
        n_pixel : int (optional, equal to 1 as default)
            Number of pixels of the output

        default : str ('disk' -> on-disk observations, 'offlimb' -> off-limb observations)

        Returns
        -------
        None
        """
        if (self.mode == 'single' and n_pixel > 1):
            raise Exception("Single pixel models cannot contain more than one pixel")

        if (default == 'disk'):
            print("Setting standard chromosphere")
            
            self.model['model'] = np.zeros((n_pixel,8), dtype=np.float64)
            self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)

            self.model['model'][:] = np.array([0.0,0.0,0.0,1.0,0.0,8.0,1.0,0.0])[None,:]
            self.model['ff'][:] = 1.0

        if (default == 'offlimb'):
            print("Setting standard chromosphere")
            
            self.model['model'] = np.zeros((n_pixel,8), dtype=np.float64)
            self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)

            self.model['model'][:] = np.array([0.0,0.0,0.0,1.0,0.0,14.0,1.0,0.0])[None,:]
            self.model['ff'][:] = 1.0
        
    def save(self, file, default=None):
        """
        Save the curent observation

        Parameters
        ----------
        file : str
            Name of the output files. Extensions will be added to it
        
        Returns
        -------
        None
        """
        
        if (self.mode == 'single'):
            print("Saving chromospheric 1D model : {0}.1d".format(file))
            f = open('{0}.1d'.format(file), 'w')            
            f.write('Bx [G]   By [G]   Bz [G]   tau    v [km/s]     deltav [km/s]   beta    a     ff\n')
            np.savetxt(f, np.atleast_2d(np.hstack([np.squeeze(self.model['model']), self.model['ff'][0]])))
            f.close()

        if (self.mode == 'multi'):
            print("Saving photospheric 3D model : {0}.h5".format(file))
            f = h5py.File('{0}.h5'.format(file), 'w')
            db_model = f.create_dataset('model', self.model['model'].shape, dtype=np.float64)
            db_ff = f.create_dataset('ff', self.model['ff'].shape, dtype=np.float64)
            db_model[:] = self.model['model']
            db_ff[:] = self.model['ff']
            f.close()