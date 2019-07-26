#from .codes import sir_code
import importlib
# from ipdb import set_trace as stop

__all__ = ["Sir"]

class Sir(object):
    
    def __init__(self):
                
        self.pysir = importlib.import_module('hazel.codes.sir_code')
        self.pysir.init()

    def synthRF(self, model, macroturbulence, filling_factor, stray):
        return self.pysir.synthRF(model, macroturbulence, filling_factor, stray)

    def synth(self, model, macroturbulence, filling_factor, stray):
        return self.pysir.synth(model, macroturbulence, filling_factor, stray)

#     def 
# def list_lines():
#     """List the lines available in SIR for synthesis
        
#     """
#     f = open('LINEAS', 'r')
#     lines = f.readlines()
#     f.close()

#     print("Available lines:")
#     for l in lines[:-1]:
#         print(l[:-1])

# def initialize(lines):
#     """Initialize the SIR synthesis code for a set of spectral lines
    
#     Args:
#         lines (list): list of lists containing the information for the lines to be synthesized
#                     Each region is defined is defined by a list containing the following elements:
#                      - A string defining which lines are synthesized in the region. E.g. '1,2,3'
#                      - Initial wavelength displacement in mA
#                      - Step in mA
#                      - Final wavelength displacement in mA

#                      E.g. lines = [['1', -500.0, 10.0, 1500.0], ['2', -750.0, 10.0, 1300.0]]
    
#     Returns:
#         int: the number of wavelength points to be synthesized

#     """
#     f = open('malla.grid', 'w')
#     f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
#     f.write("           b) The first six characters of the last line                \n")
#     f.write("          in the header (if any) must contain the symbol ---       \n")
#     f.write("\n")                                                                       
#     f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
#     f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
#     f.write("-----------------------------------------------------------------------\n")

#     for i in range(len(lines)):
#         f.write("{0}            :  {1}, {2}, {3}\n".format(lines[i][0], lines[i][1], lines[i][2], lines[i][3]))
#     f.close()

#     if (not os.path.exists('LINEAS')):
#         local = str(__file__).split('/')
#         sdir = '/'.join(local[0:-2])+'/data'
#         shutil.copy(sdir+'/LINEAS', os.getcwd())

#     if (not os.path.exists('THEVENIN')):
#         local = str(__file__).split('/')
#         sdir = '/'.join(local[0:-2])+'/data'
#         shutil.copy(sdir+'/THEVENIN', os.getcwd())
        
#     return init()

# def _interpolateNodes(logTau, nodes, nodes_logtau=None, variable=None):
#     n = logTau.shape[0]

#     if (variable is None):
#         variable = np.zeros_like(logTau)

#     out = np.zeros_like(logTau)

#     if (nodes_logtau is not None):
#         ind = np.argsort(nodes_logtau)
#         spl = PchipInterpolator(nodes_logtau[ind], nodes[ind], extrapolate=True)
#         out = spl(logTau)
#         return out

#     if (len(nodes) == 1):
#         out = variable + nodes
    
#     if (len(nodes) >= 2):
#         pos = np.linspace(0, n-1, len(nodes), dtype=int)
#         coeff = np.polyfit(logTau[pos], nodes, len(nodes)-1)
#         out = variable + np.polyval(coeff, logTau)
#     return out

# def build_model(logTau, nodes_logtau=None, nodes_T=None, nodes_vmic=None, nodes_B=None, nodes_v=None, nodes_thB=None, nodes_phiB=None,
#     var_T=None, var_vmic=None, var_B=None, var_v=None, var_thB=None, var_phiB=None):
#     """Build a SIR model given the nodes for all quantities

#     Args:
#         logTau (float): array with the output log(tau) axis
#         nodes_logtau (float): log(tau) position of the nodes
#         nodes_T (list): list with the number of nodes for Temperature [K]
#         nodes_vmic (list): list with the number of nodes for microturbulent velocity [km/s]
#         nodes_B (list): list with the number of nodes for magnetic field strength [G]
#         nodes_v (list): list with the number of nodes for velocity [km/s]
#         nodes_thB (list): list with the number of nodes for inclination of magnetic field [deg]
#         nodes_phiB (list): list with the number of nodes for azimuth of magnetic field [deg]
#         var_T (list): list with the temperature to be added to the interpolation of the nodes

#     Returns:
#         model (float): [nDepth x 7] array appropriate for synthesizing with SIR
#     """
#     n = len(logTau)
#     nodes = [nodes_T, nodes_vmic, nodes_B, nodes_v, nodes_thB, nodes_phiB]
#     variable = [var_T, var_vmic, var_B, var_v, var_thB, var_phiB]
#     model = np.zeros((n,6))
#     for i in range(6):
#         if (nodes[i] is not None):
#             if (variable[i] is None):
#                 model[:,i] = _interpolateNodes(logTau, np.asarray(nodes[i]), nodes_logtau=nodes_logtau)
#             else:
#                 model[:,i] = _interpolateNodes(logTau, np.asarray(nodes[i]), nodes_logtau=nodes_logtau, variable=variable[i])

#     return model

# def set_PSF(xPSF, yPSF):
#     """Define the spectral PSF to be convolved with the profiles
    
#     Args:
#         xPSF (float): wavelength with respect to the maximum (in mA) [nLambda]
#         yPSF (float): transmission for each wavelength displacement. It is not necessary to normalize it to unit area [nLambda]    
#     """
#     setPSF(xPSF, yPSF)

# def synthesize(model, macroturbulence=0.0, fillingFactor=1.0, stray=0.0, returnRF=True):
#     """Carry out the synthesis and returns the Stokes parameters and the response 
#     functions to all physical variables at all depths
    
#     Args:
#         model (float array): an array of size [nDepth x 7] or [nDepth x 8], where the columns contain the depth stratification of:
#             - log tau
#             - Temperature [K]
#             - Electron pressure [dyn cm^-2]  (optional)
#             - Microturbulent velocity [km/s]
#             - Magnetic field strength [G]
#             - Line-of-sight velocity [km/s]
#             - Magnetic field inclination [deg]
#             - Magnetic field azimuth [deg]
#         macroturbulence (float, optional): macroturbulence velocity [km/s]. Default: 0
#         fillingFactor (float, optional): filling factor. Default: 1
#         stray (float, optional): stray light in %. Default: 0
#         returnRF (bool, optional): return response functions
    
#     Returns:
#         stokes: (float array) Stokes parameters, with the first index containing the wavelength displacement and the remaining
#                                 containing I, Q, U and V. Size (5,nLambda)
#         rf: (float array) Response functions to T, Pe, vmic, B, v, theta, phi, all of size (4,nLambda,nDepth), plus the RF to macroturbulence of size (4,nLambda)
#                         It is not returned if returnRF=False
#     """
#     if (model.shape[1] == 7):
#         model = np.insert(model, 2, -np.ones(model.shape[0]), axis=1)
# # Boundary condition for Pe
#         model[-1,2] = 1.11634e-01        

#     if (returnRF):
#         stokes, rf = synthRF(model, macroturbulence, fillingFactor, stray)
#         return stokes, rf
#     else:
#         stokes = synth(model, macroturbulence, fillingFactor, stray)        
#         return stokes    
#     return stokes, rf

