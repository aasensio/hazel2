import numpy as np
# from ipdb import set_trace as stop

__all__ = ['hsra']

def hsra():
    """
    
    """
    path = str(__file__).split('/')
    filename = '/'.join(path[0:-1])+'/data/hsra.mod'
    f = open(filename, 'r')
    f.readline()
    ff = float(f.readline())
    f.close()
    return np.loadtxt(filename, skiprows=4), ff