import numpy as np
import pyhazel

__all__ = ["hazel"]

class hazel:
    """
    Class that synthesizes spectral lines using Hazel
    To use it:

    out = hazel()
    l, stokes, dStokes = hazel.synthDerivatives(whichPars, 1e-3, *args)
    with *args a tuple containing all parameters described in the synth method

    """

    def __init__(self):
        pyhazel._init()

    def synth(self, *args):     
        """
        Carry out a synthesis with Hazel (see the manual for the meaning of all of them)
        
        Args: 
            - synModeInput: (int, default 5) synthesis mode. Default: 5
            - nSlabsInput: (int, default 1) number of slabs. Default: 1
            - B1Input: (float, default [0,0,0]) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
            - B2Input: (float, default [0,0,0]) vector of size 3 with the magnetic field vector in spherical coordinates for the second component
            - hInput: (float, default 3.0) height
            - tau1Input: (float, default 1.0) optical depth of the first component
            - tau2Input: (float, default 0.0) optical depth of the second component        
            - boundaryInput: (float, default [0,0,0,0]) vector of size 4 with the boundary condition for (I,Q,U,V)
            - transInput: (int, default 1) transition to compute from the model atom
            - atomicPolInput: (int, default 1) include or not atomic polarization
            - magoptInput: (int, default 1) include or not magneto-optical effects
            - anglesInput: (float, default [0,0,0]) vector of size 3 describing the LOS
            - lambdaAxisInput: (float, default (-1.5,2.5) with 128 steps) vector of size 2 defining the left and right limits of the wavelength axis
            - nLambdaInput: (int, default 128) number of wavelength points
            - dopplerWidth1Input: (float, default 5) Doppler width of the first component
            - dopplerWidth2Input: (float, default 5) Doppler width of the second component
            - dampingInput: (float, default 0.0) damping
            - dopplerVelocityInput: (float, default 0.0) bulk velocity affecting the first component
            - dopplerVelocity2Input: (float, default 0.0) bulk velocity affecting the second component
            - ffInput: (float, default 1) filling factor
            - betaInput: (float, default 1) value to be multiplied by the source function of the second component to allow for emission lines in the disk
            - beta2Input: (float, default 1) enhancement factor for the source function of component 2 to allow for emission lines in the disk
            - nbarInput: (float, default [0,0,0,0]) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
            - omegaInput: (float, default [0,0,0,0]) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
            - normalization: (int, default 0) normalization of the output Stokes parameters (0-> I_max, 1-> I_peak)
		    - deltaCollision: (float, default 0.0) depolarizing rate for lower term
            
        Returns:
            - wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
            - stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
            - epsOutput: (float) array of size (4,nLambdaInput) with the emissivity vector at each wavelength
            - etaOutput: (float) array of size (4,4,nLambdaInput) with the propagation matrix at each wavelength
        """             
        return pyhazel._synth(*args)

    def __synthPerturbation(self, paramToModify, perturbation, *args):

        newPars = list(args)

        if (paramToModify == 'B1'):
            newPars[2][0] += perturbation
        if (paramToModify == 'thB1'):
            newPars[2][1] += perturbation
        if (paramToModify == 'chiB1'):
            newPars[2][2] += perturbation

        if (paramToModify == 'B2'):
            newPars[3][0] += perturbation
        if (paramToModify == 'thB2'):
            newPars[3][1] += perturbation
        if (paramToModify == 'chiB2'):
            newPars[3][2] += perturbation

        if (paramToModify == 'tau1'):
            newPars[5] += perturbation
        if (paramToModify == 'tau2'):
            newPars[6] += perturbation

        if (paramToModify == 'vth1'):
            newPars[14] += perturbation
        if (paramToModify == 'vth2'):
            newPars[15] += perturbation

        if (paramToModify == 'a'):
            newPars[16] += perturbation

        if (paramToModify == 'v1'):
            newPars[17] += perturbation
        if (paramToModify == 'v2'):
            newPars[17] += perturbation

        if (paramToModify == 'ff'):
            newPars[19] += perturbation

        if (paramToModify == 'beta1'):
            newPars[20] += perturbation
        if (paramToModify == 'beta2'):
            newPars[21] += perturbation

        newPars = tuple(newPars)

        _, stokesPerturbed, _, _ = pyhazel._synth(*newPars)

        return stokesPerturbed


    def synth_RF(self, paramsRF, perturbation=1e-3, *args):
        """
        Carry out a synthesis with Hazel and also compute response functions to the parameters
        
        Args: 
            (see the manual for the meaning of all of them)
            - paramsRF : (str) list containing any combination of the following parameters to which the RF are computed
                ['B1', 'thB1', 'chiB1', 'B2', 'thB2', 'chiB2', 'tau1', 'tau2', 'vth1','vth2','a','v1','v2', 'ff', 'beta', 'beta2']
            - perturbation: (float, default 1e-3) perturbation used for computing the numerical RF
            - synModeInput: (int, default 5) synthesis mode. Default: 5
            - nSlabsInput: (int, default 1) number of slabs. Default: 1
            - B1Input: (float, default [0,0,0]) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
            - B2Input: (float, default [0,0,0]) vector of size 3 with the magnetic field vector in spherical coordinates for the second component
            - hInput: (float, default 3.0) height
            - tau1Input: (float, default 1.0) optical depth of the first component
            - tau2Input: (float, default 0.0) optical depth of the second component        
            - boundaryInput: (float, default [0,0,0,0]) vector of size 4 with the boundary condition for (I,Q,U,V)
            - transInput: (int, default 1) transition to compute from the model atom
            - atomicPolInput: (int, default 1) include or not atomic polarization
            - magoptInput: (int, default 1) include or not magneto-optical effects
            - anglesInput: (float, default [0,0,0]) vector of size 3 describing the LOS
            - lambdaAxisInput: (float, default (-1.5,2.5) with 128 steps) vector of size 2 defining the left and right limits of the wavelength axis
            - nLambdaInput: (int, default 128) number of wavelength points
            - dopplerWidth1Input: (float, default 5) Doppler width of the first component
            - dopplerWidth2Input: (float, default 5) Doppler width of the second component
            - dampingInput: (float, default 0.0) damping
            - dopplerVelocityInput: (float, default 0.0) bulk velocity affecting the first component
            - dopplerVelocity2Input: (float, default 0.0) bulk velocity affecting the second component
            - ffInput: (float, default 1) filling factor
            - betaInput: (float, default 1) value to be multiplied by the source function of the second component to allow for emission lines in the disk
            - beta2Input: (float, default 1) enhancement factor for the source function of component 2 to allow for emission lines in the disk
            - nbarInput: (float, default [0,0,0,0]) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
            - omegaInput: (float, default [0,0,0,0]) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
            - normalization: (int, default 0) normalization of the output Stokes parameters (0-> I_max, 1-> I_peak)
		    - deltaCollision: (float, default 0.0) depolarizing rate for lower term
        Returns:
            - wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
            - stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
            - stokesDeriv: (float) array of size (nPar,4,nLambdaInput) with the response functions or each indicated parameter
        """ 

        wavelengthOutput, stokes, epsOutput, etaOutput = pyhazel._synth(*args)

        nRF = len(paramsRF)

        stokesDeriv = np.zeros((nRF,4,len(wavelengthOutput)))

        for i in range(nRF):
            stokesNew = self.__synthPerturbation(paramsRF[i], perturbation, *args)

            stokesDeriv[i,:,:] = (stokesNew - stokes) / perturbation            

        return wavelengthOutput, stokes, stokesDeriv
