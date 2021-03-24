import numpy as np
import h5py

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
        
        self.model = {'model': None, 'ff': None, 'vmac': None}

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
        self.model['vmac'] = np.zeros((n_pixel,), dtype=np.float64)

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
        
        print("Setting photosphere to model {0}".format(default))
        path = str(__file__).split('/')
        filename = '/'.join(path[0:-1])+'/data/{0}.1d'.format(default)
        f = open(filename, 'r')
        f.readline()
        tmp = f.readline().split()
        ff, vmac = float(tmp[0]), float(tmp[1])
        f.close()
        model = np.loadtxt(filename, skiprows=4)

        nz = model.shape[0]

        self.model['model'] = np.zeros((n_pixel,nz,8), dtype=np.float64)
        self.model['ff'] = np.zeros((n_pixel,), dtype=np.float64)
        self.model['vmac'] = np.zeros((n_pixel,), dtype=np.float64)

        self.model['model'][:] = model[None,:,:]
        self.model['ff'][:] = ff
        self.model['vmac'][:] = vmac

    def list_models(self):
        docs = """
        All models have been extracted from SIR (https://github.com/BasilioRuiz/SIR-code/tree/master/models)
        
        One-component quiet Sun models:

            holmu11.mod ...  Holweger H., & Muller E.A., 1974, Sol Phys. 39 19

            hsra11.mod   ...  Harvard Smithsonian Reference Atmosphere (Gingerich O.,
                            Noyes R.W., Kalkofen W., & Cuny Y., 1971. Sol. Phys. 18, 347)

            valc11.mod   ...  Vernaza J.E., Avrett E.H., & Loeser R., 1981, ApJS 45, 635

            mackkl11.mod ...  Maltby P., Avrett E.H., Carlsson M., Kjeldseth-Moe O.,
                            Kurucz R.L., & Loeser R., 1986 ApJ 306 284

            nelsoncold.mod  Nelson G. P. 1978 Sol. Phys. 60, 5

            nelsonhot.mod   Nelson G. P. 1978 Sol. Phys. 60, 5

            grevesse11.mod .. Grevesse N., Sauval A.J. 1999 A&A 347, 348
            

        Sunspot Models:

            emaltby11.mod .. (E model) Maltby P., Avrett E.H., Carlsson M., 
                            Kjeldseth-Moe O., Kurucz R.L., & Loeser R., 1986 ApJ 306 284
                            
            mmaltby11.mod ..  (M model) Maltby P., Avrett E.H., Carlsson M., 
                            Kjeldseth-Moe O., Kurucz R.L., & Loeser R., 1986 ApJ 306 284
                            

            cool11.mod ...    Collados M., Martínez Pillet V., Ruiz Cobo B., 
                            Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 
                            (Umbral model for a big spot)

            hot11.mod  ...    Collados M., Martínez Pillet V., Ruiz Cobo B., 
                            Del Toro Iniesta J.C., & Vázquez M. 1994 A&A 291 622 
                    (Umbral model for a small spot)

            penumjti11.mod .. Del Toro Iniesta J.C., Tarbell T.D., Ruiz Cobo B., 
                            1994, ApJ 436,400 
                    (penumbral model)

        Active Regions Models:

            solannt11.mod ..  (network model, S. K. Solanki, private comunication)

            solanpl11.mod ..  (plage model, S. K. Solanki, private comunication)

            Fontenla Models:

        Photospheric models
            FALB = MODEL1001

            FALC11 - QS model cell center
            An area with the same intensity as the median in a histogram of a Ca II K image of a quiet area of the Sun. 
            We find that the median is very   close to the peak of the distribution but is statistically more stable. 
            These intensities correspond to most of the central area of supergranular cells that are usually known as "quiet-Sun cell interior."

            FALD = MODEL1002  : Network
            
            FALE11- QS model network
            A bright area separating supergranulation (or network) cells, often called "network lane." We describe this category as "quiet-Sun network."

            FALF = MODEL1003 - QS model active network
                Certain network lane areas that are much brighter than the average. We describe this category as "active network."
                It spans from logtau=1.3 to -6.8

            FALF11=  FALF but spans from logtau=1.4 to -5.7
            
            FALH = MODEL1004-  Plage model
                It spans from logtau=1.5 to -6.6
                
            FALH11= FALH but spans from logtau=1.4 to -5.7
            
            FALP= MODEL1005- Facula model
                    It spans from logtau=1.3 to -6.7
                    
            FALP11=FALP but spans from logtau=1.1 to -5.7
            
            FALR11 - penumbral model
            
            FALS= MODEL1006 - sunspot model spans from logtau=1.4 to -5.4
            
            FALS11= FALS  but spans from logtau=1.3 to -4.9
        """
        print(docs)
        
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
            f.write('ff  vmac\n')
            f.write('{0}  {1}\n'.format(self.model['ff'][0], self.model['vmac'][0]))
            f.write('\n')
            f.write('  logtau     T        Pe           vmic        v            Bx           By         Bz\n')
            
            np.savetxt(f, np.squeeze(self.model['model']))
            f.close()

        if (self.mode == 'multi'):
            print("Saving photospheric 3D model : {0}.h5".format(file))
            f = h5py.File('{0}.h5'.format(file), 'w')
            db_model = f.create_dataset('model', self.model['model'].shape, dtype=np.float64)
            db_ff = f.create_dataset('ff', self.model['ff'].shape, dtype=np.float64)
            db_vmac = f.create_dataset('vmac', self.model['vmac'].shape, dtype=np.float64)
            db_model[:] = self.model['model']
            db_ff[:] = self.model['ff']
            db_vmac[:] = self.model['vmac']
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