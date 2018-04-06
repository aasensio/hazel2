import numpy as np
import h5py
from astropy.io import fits
import os
# from ipdb import set_trace as stop

__all__ = ['Generic_output_file', 'Generic_observed_file', 'Generic_hazel_file', 'Generic_SIR_file', 'Generic_parametric_file', 'Generic_stray_file']

class Generic_output_file(object):    

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

        if (self.extension == '1d'):
            raise Exception("1d files not allowed as output")

    def open(self, model):        
        if (self.extension == 'h5'):
            # Open the file            
            self.handler = h5py.File(self.filename, 'w')
            
            # Generate all handlers for things we'll write here
            self.out_spectrum = {}
            for k, v in model.spectrum.items():
                n_stokes, n_lambda = v.stokes.shape
                self.out_spectrum[k] = self.handler.create_dataset(k, (model.n_pixels, n_stokes, n_lambda))

            if (model.working_mode == 'inversion'):
                self.out_model = {}
                for k, v in model.atmospheres.items():
                    db = self.handler.create_group(k)
                    self.out_model[k] = {}
                    for k2, v2 in v.reference.items():                    
                        if (np.isscalar(v2)):
                            n_depth = 1
                        else:
                            n_depth = len(v2)
                        self.out_model[k][k2] = db.create_dataset(k2, (model.n_pixels, n_depth))
            
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def write(self, model, pixel=None):

        if (self.extension == 'h5'):
            for k, v in model.spectrum.items():                
                self.out_spectrum[k][pixel,...] = v.stokes

            if (model.working_mode == 'inversion'):
                for k, v in model.atmospheres.items():                
                    for k2, v2 in v.reference.items():                    
                        self.out_model[k][k2][pixel,...] = v2

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

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            f = open(self.filename, 'r')
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

        if (self.extension == 'h5'):
            los = self.handler['LOS'][pixel,...]
            mu = np.cos(los[0] * np.pi / 180.0)
            return self.handler['stokes'][pixel,...].T, self.handler['sigma'][pixel,...].T, los, mu, self.handler['boundary'][pixel,...].T

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
            tmp, _, _ = self.handler['stokes'].shape
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

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            tmp = np.loadtxt(self.filename, skiprows=1)
            return tmp

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

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

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            return np.loadtxt(self.filename, skiprows=1)

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

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

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            return np.loadtxt(self.filename, skiprows=1)

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...], self.handler['ff'][pixel]

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
        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            f = open(self.filename, 'r')
            f.readline()
            ff = float(f.readline())
            f.close()
            return np.loadtxt(self.filename, skiprows=4), ff

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...], self.handler['ff'][pixel]

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
            tmp, _, _ = self.handler['model'].shape
            self.close()
            return tmp