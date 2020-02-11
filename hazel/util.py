import numpy as np
import h5py
from asciitree import LeftAligned
from collections import OrderedDict
from asciitree.drawing import BoxStyle, BOX_DOUBLE, BOX_BLANK
# from ipdb import set_trace as stop

__all__ = ['i0_allen', '_extract_parameter_cycles', 'isint', 'fvoigt', 'lower_dict_keys', 'show_tree']

def i0_allen(wavelength, muAngle):
    """
    Return the solar intensity at a specific wavelength and heliocentric angle
    wavelength: wavelength in angstrom
    muAngle: cosine of the heliocentric angle
    """
    C = 2.99792458e10
    H = 6.62606876e-27

    if (muAngle == 0):
        return 0.0

    lambdaIC = 1e4 * np.asarray([0.20,0.22,0.245,0.265,0.28,0.30,0.32,0.35,0.37,0.38,0.40,0.45,0.50,0.55,0.60,0.80,1.0,1.5,2.0,3.0,5.0,10.0])
    uData = np.asarray([0.12,-1.3,-0.1,-0.1,0.38,0.74,0.88,0.98,1.03,0.92,0.91,0.99,0.97,0.93,0.88,0.73,0.64,0.57,0.48,0.35,0.22,0.15])
    vData = np.asarray([0.33,1.6,0.85,0.90,0.57, 0.20, 0.03,-0.1,-0.16,-0.05,-0.05,-0.17,-0.22,-0.23,-0.23,-0.22,-0.20,-0.21,-0.18,-0.12,-0.07,-0.07])

    lambdaI0 = 1e4 * np.asarray([0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.48,0.50,0.55,0.60,0.65,0.70,0.75,\
        0.80,0.90,1.00,1.10,1.20,1.40,1.60,1.80,2.00,2.50,3.00,4.00,5.00,6.00,8.00,10.0,12.0])
    I0 = np.asarray([0.06,0.21,0.29,0.60,1.30,2.45,3.25,3.77,4.13,4.23,4.63,4.95,5.15,5.26,5.28,5.24,5.19,5.10,5.00,4.79,4.55,4.02,3.52,3.06,2.69,2.28,2.03,\
        1.57,1.26,1.01,0.81,0.53,0.36,0.238,0.160,0.078,0.041,0.0142,0.0062,0.0032,0.00095,0.00035,0.00018])
    I0 *= 1e14 * (lambdaI0 * 1e-8)**2 / C

    u = np.interp(wavelength, lambdaIC, uData)
    v = np.interp(wavelength, lambdaIC, vData)
    i0 = np.interp(wavelength, lambdaI0, I0)
    
    return (1.0 - u - v + u * muAngle + v * muAngle**2)* i0

def _extract_parameter_cycles(s):
    tmp = s[0].split('->')
    value = float(tmp[0])
    cycle1 = tmp[1].strip()
    cycles = [cycle1] + s[1:]

    return value, cycles

def isint(str):
    try:
        int(str)
        return True
    except ValueError:
        return False

def isfloat(str):
    if (str is None):
        return False
    try:
        float(str)
        return True
    except ValueError:
        return False

def toint(l):
    return [int(x) if isint(x) else x for x in l]

def tofloat(l):
    return [float(x) if isfloat(x) else None for x in l]

def tobool(l):
    return True if l == 'True' else False

def onlyint(l):
    return [i for i in l if isinstance(i, int)]


def fvoigt(damp,v):
    
    """
    Fast implementation of the Voigt-Faraday function

    Parameters
    ----------
        damp : float
            damping parameter

        v : float
            normalized wavelength (lambda-lambda0) / sigma
        
    Returns
    -------
        voigt, faraday : float
            Value of the Voigt and Faraday functions


    Notes
    ----- 
        A rational approximation to the complex error function is used
        after Hui, Armstrong, and Wray(1978, JQSRT 19, 509). H and F are 
        the real and imaginary parts of such function, respectively.
        The procedure is inspired on that in SIR (Ruiz Cobo & del Toro 
        Iniesta 1992, ApJ 398, 385). On its turn, that routine was taken
        from modifications by A. Wittmann (1986) to modifications by S.K.
        Solanki (1985) to an original FORTRAN routine written by J.W. Harvey
        and A. Nordlund.
    """
    
    A = [122.607931777104326, 214.382388694706425, 181.928533092181549,\
        93.155580458138441, 30.180142196210589, 5.912626209773153,\
        0.564189583562615]

    B = [122.60793177387535, 352.730625110963558, 457.334478783897737,\
        348.703917719495792, 170.354001821091472, 53.992906912940207,\
        10.479857114260399,1.]

    z = np.array(damp*np.ones(len(v)) + -abs(v)*1j)

    Z = ((((((A[6]*z+A[5])*z+A[4])*z+A[3])*z+A[2])*z+A[1])*z+A[0])/\
    (((((((z+B[6])*z+B[5])*z+B[4])*z+B[3])*z+B[2])*z+B[1])*z+B[0])

    h = np.real(Z)
    f = np.sign(v)*np.imag(Z)*0.5

    return h, f


def lower_dict_keys(d):
    out = {}
    for k, v in d.items():
        out[k.lower()] = v
    return out

def show_tree(hdf5_file):
    tree = {hdf5_file: OrderedDict()}

    f = h5py.File(hdf5_file, 'r')
    for k, v in f.items():
        tree[hdf5_file][k] = OrderedDict()
        for k2, v2 in v.items():
            tree[hdf5_file][k][f'{k2} -> {v2.shape}  {v2.dtype}'] = OrderedDict()            

    chrs = dict(
            UP_AND_RIGHT=u"\u2514",
            HORIZONTAL=u"\u2500",
            VERTICAL=u"\u2502",
            VERTICAL_AND_RIGHT=u"\u251C"
        )

    tr = LeftAligned(draw=BoxStyle(gfx = chrs, horiz_len=1))
    print(tr(tree))