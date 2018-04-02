Python wrapper
==============


We have developed a wrapper so that the code can be called from Python easily
and synthesize Stokes profiles using a simple interface.
The wrapper is installed by going to the ``runPy`` directory and 
typing 

::

    python setup.py build_ext --inplace


It will generate a ``pyhazel.so`` file appropriate for your architecture, that
can be imported from Python. The ``test.py`` file shows how to call the
wrapper to synthesize the Stokes profiles. In the following we describe the
inputs:

::

    Args: (see the manual for the meaning of all of them)
        synModeInput: (int) synthesis mode
        nSlabsInput: (int) number of slabs
        B1Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
        B2Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the second component
        hInput: (float) height
        tau1Input: (float) optical depth of the first component
        tau2Input: (float) optical depth of the second component        
        boundaryInput: (float) vector of size 4 with the boundary condition for (I,Q,U,V)
        transInput: (int) transition to compute from the model atom
        atomicPolInput: (int) include or not atomic polarization
        anglesInput: (float) vector of size 3 describing the LOS
        lambdaAxisInput: (float) vector of size 2 defining the left and right limits of the wavelength axis
        nLambdaInput: (int) number of wavelength points
        dopplerWidth1Input: (float) Doppler width of the first component
        dopplerWidth2Input: (float) Doppler width of the second component
        dampingInput: (float) damping
        dopplerVelocityInput: (float) bulk velocity affecting the first component
        dopplerVelocity2Input: (float) bulk velocity affecting the second component
        ffInput: (float) filling factor
        betaInput: (float) enhancement factor for the source function of component 1 to allow for emission lines in the disk
        beta2Input: (float) enhancement factor for the source function of component 2 to allow for emission lines in the disk
        nbarInput: (float) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
        omegaInput: (float) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
        normalization: (int) normalization of the output Stokes parameters (0-> I_max, 1-> I_peak)
        deltaCollision: (float) depolarizing rate for lower term
        
    Returns:
        wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
        stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
        epsOutput: (float) array of size (4,nLambdaInput) with the emissivity vector at each wavelength
        etaOutput: (float) array of size (4,4,nLambdaInput) with the propagation matrix at each wavelength
    """
