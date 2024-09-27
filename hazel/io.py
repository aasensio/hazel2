import numpy as np
import h5py
from astropy.io import fits
import hazel
import os
import datetime
import warnings

try:
    import zarr
except:
    pass
    # warnings.warn("zarr module not found. You will not be able to use zarr as input/output.")
    #EDGAR: add your ncdf routines here

__all__ = ['Generic_output_file', 'Generic_observed_file', 'Generic_hazel_file', 'Generic_SIR_file', 'Generic_parametric_file', 'Generic_stray_file', 'Generic_mask_file']

class Generic_output_file(object):    

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

        if (self.extension == '1d'):
            raise Exception("1d files not allowed as output")

    def open(self, model):

        root = model.root        
        
        if (self.extension == 'h5' or self.extension == 'zarr'):
            
            # Open the file        
            if (self.extension == 'h5'):                
                self.handler = h5py.File(self.filename, 'w')
                
            if (self.extension == 'zarr'):
                self.handler = zarr.open(self.filename, 'w')

            self.handler.attrs['version'] = hazel.__version__
            self.handler.attrs['date'] = datetime.datetime.today().isoformat(' ')

            # Save configuration file
            self.handler.attrs['configuration'] = [a.encode('utf8') for a in model.configuration.configuration_txt]
            
            # Generate all handlers for things we'll write here
            self.out_spectrum = {}
            for k, v in model.spectrum.items():
                db = self.handler.create_group(k)                
                
                self.out_spectrum[k] = {}

                n_stokes, n_lambda = v.stokes.shape
                self.out_spectrum[k]['wavelength'] = db.create_dataset('wavelength', shape=(n_lambda,), dtype=np.float64)                
                                
                if (model.save_all_cycles):
                    self.out_spectrum[k]['stokes'] = db.create_dataset('stokes', shape=(model.n_pixels, model.n_randomization, model.n_cycles, n_stokes, n_lambda), dtype=np.float64)
                    self.out_spectrum[k]['stokes'].dims[0].label = 'pixel'
                    self.out_spectrum[k]['stokes'].dims[1].label = 'randomization'
                    self.out_spectrum[k]['stokes'].dims[2].label = 'cycle'
                    self.out_spectrum[k]['stokes'].dims[3].label = 'stokes_parameter'
                    self.out_spectrum[k]['stokes'].dims[4].label = 'wavelength'
                else:
                    self.out_spectrum[k]['stokes'] = db.create_dataset('stokes', shape=(model.n_pixels, model.n_randomization, n_stokes, n_lambda), dtype=np.float64)
                    self.out_spectrum[k]['stokes'].dims[0].label = 'pixel'
                    self.out_spectrum[k]['stokes'].dims[1].label = 'randomization'                    
                    self.out_spectrum[k]['stokes'].dims[2].label = 'stokes_parameter'
                    self.out_spectrum[k]['stokes'].dims[3].label = 'wavelength'

                if (v.interpolate_to_lr):
                    n_stokes, n_lambda = v.stokes_lr.shape
                    
                    self.out_spectrum[k]['wavelength_lr'] = db.create_dataset('wavelength_lr', shape=(n_lambda,), dtype=np.float64)                

                    if (model.save_all_cycles):
                        self.out_spectrum[k]['stokes_lr'] = db.create_dataset('stokes_lr', shape=(model.n_pixels, model.n_randomization, model.n_cycles, n_stokes, n_lambda), dtype=np.float64)
                        self.out_spectrum[k]['stokes_lr'].dims[0].label = 'pixel'
                        self.out_spectrum[k]['stokes_lr'].dims[1].label = 'randomization'
                        self.out_spectrum[k]['stokes_lr'].dims[2].label = 'cycle'
                        self.out_spectrum[k]['stokes_lr'].dims[3].label = 'stokes_parameter'
                        self.out_spectrum[k]['stokes_lr'].dims[4].label = 'wavelength'
                    else:
                        self.out_spectrum[k]['stokes_lr'] = db.create_dataset('stokes_lr', shape=(model.n_pixels, model.n_randomization, n_stokes, n_lambda), dtype=np.float64)
                        self.out_spectrum[k]['stokes_lr'].dims[0].label = 'pixel'
                        self.out_spectrum[k]['stokes_lr'].dims[1].label = 'randomization'                    
                        self.out_spectrum[k]['stokes_lr'].dims[2].label = 'stokes_parameter'
                        self.out_spectrum[k]['stokes_lr'].dims[3].label = 'wavelength'

                if (model.working_mode == 'inversion'):

                    self.out_spectrum[k]['chi2'] = db.create_dataset('chi2', shape=(model.n_pixels, model.n_randomization, model.n_cycles), dtype=np.float64)
                    self.out_spectrum[k]['chi2'].dims[0].label = 'pixel'
                    self.out_spectrum[k]['chi2'].dims[1].label = 'randomization'
                    self.out_spectrum[k]['chi2'].dims[2].label = 'cycle'
                
                    self.out_spectrum[k]['bic'] = db.create_dataset('bic', shape=(model.n_pixels, model.n_randomization, model.n_cycles), dtype=np.float64)
                    self.out_spectrum[k]['bic'].dims[0].label = 'pixel'
                    self.out_spectrum[k]['bic'].dims[1].label = 'randomization'
                    self.out_spectrum[k]['bic'].dims[2].label = 'cycle'

                    self.out_spectrum[k]['aic'] = db.create_dataset('aic', shape=(model.n_pixels, model.n_randomization, model.n_cycles), dtype=np.float64)
                    self.out_spectrum[k]['aic'].dims[0].label = 'pixel'
                    self.out_spectrum[k]['aic'].dims[1].label = 'randomization'
                    self.out_spectrum[k]['aic'].dims[2].label = 'cycle'

            if (model.working_mode == 'inversion'):
                self.out_model = {}
                self.out_error = {}
                self.out_nodes = {}
                for k, v in model.atmospheres.items():
                    db = self.handler.create_group(k)
                    self.out_model[k] = {}
                    self.out_error[k] = {}
                    self.out_nodes[k] = {}

                    # Save metadata for type of reference for magnetic field components
                    if (hasattr(v, 'reference_frame') or hasattr(v, 'ref frame')):
                        db.attrs['reference frame'] = v.reference_frame

                    db.attrs['ranges'] = str(dict(v.ranges))
                    db.attrs['regularization'] = str(dict(v.regularization))
                    db.attrs['cycles'] = str(dict(v.cycles))
                    
                    if (hasattr(v, 'log_tau')):
                        self.out_model[k]['log_tau'] = db.create_dataset('log_tau', shape=(len(v.log_tau),), dtype=np.float64)
                    
                    for k2, v2 in v.reference.items():
                        if (np.isscalar(v2)):
                            n_depth = 1
                        else:
                            n_depth = len(v2)                            
                        
                        if (model.save_all_cycles):
                            self.out_model[k][k2] = db.create_dataset(k2, shape=(model.n_pixels, model.n_randomization, model.n_cycles, n_depth), dtype=np.float64)
                            self.out_model[k][k2].dims[0].label = 'pixel'
                            self.out_model[k][k2].dims[1].label = 'randomization'
                            self.out_model[k][k2].dims[2].label = 'cycle'
                            self.out_model[k][k2].dims[3].label = 'depth'

                            d_nodes = h5py.special_dtype(vlen=np.dtype('float64'))
                            self.out_error[k][k2] = db.create_dataset('{0}_err'.format(k2), shape=(model.n_pixels, model.n_randomization, model.n_cycles,), dtype=d_nodes)
                            self.out_error[k][k2].dims[0].label = 'pixel'
                            self.out_error[k][k2].dims[1].label = 'randomization'
                            self.out_error[k][k2].dims[2].label = 'cycle'                        
                            
                            d_nodes = h5py.special_dtype(vlen=np.dtype('float64'))
                            self.out_nodes[k][k2] = db.create_dataset('{0}_nodes'.format(k2), shape=(model.n_pixels, model.n_randomization, model.n_cycles,), dtype=d_nodes)
                            self.out_nodes[k][k2].dims[0].label = 'pixel'
                            self.out_nodes[k][k2].dims[1].label = 'randomization'
                            self.out_nodes[k][k2].dims[2].label = 'cycle'
                        else:
                            self.out_model[k][k2] = db.create_dataset(k2, shape=(model.n_pixels, model.n_randomization, n_depth), dtype=np.float64)
                            self.out_model[k][k2].dims[0].label = 'pixel'
                            self.out_model[k][k2].dims[1].label = 'randomization'                            
                            self.out_model[k][k2].dims[2].label = 'depth'

                            d_nodes = h5py.special_dtype(vlen=np.dtype('float64'))
                            self.out_error[k][k2] = db.create_dataset('{0}_err'.format(k2), shape=(model.n_pixels, model.n_randomization,), dtype=d_nodes)
                            self.out_error[k][k2].dims[0].label = 'pixel'
                            self.out_error[k][k2].dims[1].label = 'randomization'                            
                            
                            d_nodes = h5py.special_dtype(vlen=np.dtype('float64'))
                            self.out_nodes[k][k2] = db.create_dataset('{0}_nodes'.format(k2), shape=(model.n_pixels, model.n_randomization,), dtype=d_nodes)
                            self.out_nodes[k][k2].dims[0].label = 'pixel'
                            self.out_nodes[k][k2].dims[1].label = 'randomization'
                            
                    for k2, v2 in v.units.items():                        
                        self.out_model[k][k2].attrs['unit'] = v2
                        
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def write(self, model, pixel=0, randomization=0):

        if (self.extension == 'h5'):            
            for k, v in model.spectrum.items():                            
                for cycle in range(model.n_cycles):
                    self.out_spectrum[k]['wavelength'][:] = v.wavelength_axis
                    if (model.save_all_cycles):
                        self.out_spectrum[k]['stokes'][pixel,randomization,cycle,...] = v.stokes_cycle[cycle]
                    else:
                        self.out_spectrum[k]['stokes'][pixel,randomization,...] = v.stokes_cycle[cycle]

                    if (v.interpolate_to_lr):
                        
                        self.out_spectrum[k]['wavelength_lr'][:] = v.wavelength_axis_lr

                        if (model.save_all_cycles):
                            self.out_spectrum[k]['stokes_lr'][pixel,randomization,cycle,...] = v.stokes_lr_cycle[cycle]
                        else:
                            self.out_spectrum[k]['stokes_lr'][pixel,randomization,...] = v.stokes_lr_cycle[cycle]


                    if (model.working_mode == 'inversion'):
                        if (model.save_all_cycles):
                            self.out_spectrum[k]['chi2'][pixel,randomization,cycle] = v.chi2_cycle[cycle]
                            self.out_spectrum[k]['bic'][pixel,randomization,cycle] = v.bic_cycle[cycle]
                            self.out_spectrum[k]['aic'][pixel,randomization,cycle] = v.aic_cycle[cycle]
                        else:
                            self.out_spectrum[k]['chi2'][pixel,randomization] = v.chi2_cycle[cycle]
                            self.out_spectrum[k]['bic'][pixel,randomization] = v.bic_cycle[cycle]
                            self.out_spectrum[k]['aic'][pixel,randomization] = v.aic_cycle[cycle]

            if (model.working_mode == 'inversion'):                
                for k, v in model.atmospheres.items():
                    if (hasattr(v, 'log_tau')):
                        self.out_model[k]['log_tau'][:] = v.log_tau
                    for cycle in range(model.n_cycles):

                        # Model parameters
                        for k2, v2 in v.reference_cycle[cycle].items():
                            if (model.save_all_cycles):
                                self.out_model[k][k2][pixel,randomization,cycle,...] = v2
                            else:
                                self.out_model[k][k2][pixel,randomization,...] = v2
                            
                        # Model node positions                                                                  
                        for k2, v2 in v.nodes_location_cycle[cycle].items():
                            if (model.save_all_cycles):
                                self.out_nodes[k][k2][pixel,randomization,cycle] = np.atleast_1d(v2)
                            else:
                                self.out_nodes[k][k2][pixel,randomization] = np.atleast_1d(v2)

                        # Model parameter errors
                        for k2, v2 in v.error_cycle[cycle].items():
                            if (model.save_all_cycles):
                                self.out_error[k][k2][pixel,randomization,cycle] = np.atleast_1d(v2)
                            else:
                                self.out_error[k][k2][pixel,randomization] = np.atleast_1d(v2)
                            

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):            
            self.handler.close()
            del self.handler
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp


class Generic_observed_file(object):

    def __init__(self, filename, root):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename
        self.root = root

    def open(self):
        
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.root+self.filename, 'r')

            # Check if we share sigma and boundary conditions for all pixels
            self.ndim_sigma = self.handler['sigma'].ndim
            self.ndim_boundary = self.handler['boundary'].ndim
            return

        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):            
            f = open(self.root+self.filename, 'r')
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

            return stokes, noise, los, mu, boundary[:,None]*np.ones((4,n_lambda))

        if (self.extension == 'h5' or self.extension == 'zarr'):
            los = self.handler['LOS'][pixel,...]
            mu = np.cos(los[0] * np.pi / 180.0)

            # Check if we share sigma and boundary conditions for all pixels
            if (self.ndim_sigma == 3):
                sigma = self.handler['sigma'][pixel,...].T
            else:
                sigma = self.handler['sigma'][:].T

            # Check if we share sigma and boundary conditions for all pixels
            if (self.ndim_sigma == 3):
                boundary = self.handler['boundary'][pixel,...].T
            else:
                boundary = self.handler['boundary'][:].T

            return self.handler['stokes'][pixel,...].T, sigma, los, mu, boundary

        
        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp, _, _ = self.handler['stokes'].shape
            self.close()
            return tmp

class Generic_mask_file(object):

    def __init__(self, filename, root):
        self.filename = filename
        self.root = root
        if (filename is not None):
            self.extension = os.path.splitext(filename)[1][1:]        

    def open(self):        

        if (self.filename is None):
            return

        if (self.extension == '1d'):
            raise Exception("1D files are not allowed for masks.")
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.root + self.filename, 'r')
            return

        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.root + self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.root + self.filename, memmap=True)
            return

    def read(self, pixel=None):        
        if (self.filename is None):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):            
            return self.handler['mask'][pixel]
        

    def close(self):        
        if (self.filename is None):
            return
            
        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler
        

    def get_npixel(self):        
        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp = len(self.handler['mask'])
            self.close()
            return tmp        

class Generic_stray_file(object):

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            f = open(self.filename, 'r')
            f.readline()
            v, ff = np.array(f.readline().split()).astype('float64')
            f.readline()
            f.readline()
                        
            stray_profile = np.loadtxt(self.filename, skiprows=4)
            return [stray_profile, v], ff

        if (self.extension == 'h5' or self.extension == 'zarr'):
            return [self.handler['profile'][pixel,...], self.handler['model'][pixel,...]], self.handler['ff'][pixel]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_parametric_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            tmp = np.loadtxt(self.filename, skiprows=1)
            return tmp[0:-1], tmp[-1]

        if (self.extension == 'h5' or self.extension == 'zarr'):
            return self.handler['model'][pixel,...], self.handler['ff'][pixel]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_hazel_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):            
            tmp = np.loadtxt(self.filename, skiprows=1)
            return tmp[0:-1], tmp[-1]

        if (self.extension == 'h5' or self.extension == 'zarr'):            
            return self.handler['model'][pixel,...], self.handler['ff'][pixel]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_SIR_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return
        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return
        if (self.extension == 'zarr'):
            self.handler = zarr.open(self.filename, 'r')
            return
        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            f = open(self.filename, 'r')
            f.readline()
            tmp = f.readline().split()
            ff, vmac = float(tmp[0]), float(tmp[1])
            f.close()
            return np.loadtxt(self.filename, skiprows=4), ff, vmac

        if (self.extension == 'h5' or self.extension == 'zarr'):            
            return self.handler['model'][pixel,...], self.handler['ff'][pixel], self.handler['vmac'][pixel]

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.handler.close()
            del self.handler

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5' or self.extension == 'zarr'):
            self.open()
            tmp, _, _ = self.handler['model'].shape
            self.close()
            return tmp
