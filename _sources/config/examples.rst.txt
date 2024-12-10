.. include:: ../hazel_name
.. _configuration:

Examples
========

There is a single human-readable configuration file that controls the behavior of the code. It can be
used both for synthesis and inversion.


Synthesis
---------

In the following we paste a typical configuration file that can be used for synthesis. Later we describe all sections.
As you can see, you first need to define the working mode, then the spectral regions of interest and finally the
specific details for each atmosphere.

::

    # |hazel2| configuration File

    [Working mode]
    Output file = output.h5
    

    [Spectral regions]
        [[Region 1]]
        Name = spec1
        Wavelength = 10826, 10833, 150
        Topology = ph2 -> ch1 -> te1 -> st1    
        LOS = 0.0, 0.0, 90.0
        Boundary condition = 1.0, 0.0, 0.0, 0.0       # I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)
        Wavelength file = 'observations/10830.wavelength'
        Wavelength weight file = 'observations/10830.weights'
        Instrumental profile = 3.0
    
    [Atmospheres]

        [[Photosphere 2]]
        Name = ph2
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300,
        Reference frame = vertical
    
        [[Chromosphere 1]]
        Name = ch1                                              # Name of the atmosphere component
        Spectral region = spec1                                 # Spectral region to be used for synthesis
        Height = 3.0                                            # Height of the slab
        Line = 10830                                            # 10830, 5876
        Wavelength = 10826, 10833                         # Wavelength range used for synthesis
        Reference atmospheric model = 'chromospheres/model_chromosphere.1d'    # File with model parameters
        Reference frame = vertical

        [[Parametric 1]]
        Name = te1
        Spectral region = spec1
        Wavelength = 10826, 10833
        Reference atmospheric model = 'telluric/model_telluric.1d'    # File with model parameters
        Type = Gaussian           # Gaussian, Voigt, MoGaussian, MoVoigt 

        [[Straylight 1]]
        Name = st1
        Spectral region = spec1
        Wavelength = 10826, 10833    
        Reference atmospheric model = 'straylight/model_stray.1d'    # File with model parameters

Inversion
---------

In the following we paste a typical configuration file for the inversion mode.

::

    # |hazel2| configuration File

    [Working mode]
    Output file = output.h5
    Number of cycles = 1
    Maximum iterations = 10
    Relative error = 1e-4

    [Spectral regions]
        [[Region 1]]
        Name = spec1
        Wavelength = 10826, 10833, 150
        Topology = ph2 -> ch1 -> te1 #-> st1    
        LOS = 0.0, 0.0, 90.0
        Boundary condition = 1.0, 0.0, 0.0, 0.0       # I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)
        Wavelength file = 'observations/10830.wavelength'
        Wavelength weight file = 'observations/10830.weights'
        Observations file = 'observations/10830_stokes.h5'
        Straylight file = 'observations/10830_stray.1d'
        Mask file = None
        Weights Stokes I = 1.0, 1.0, 1.0, 1.0
        Weights Stokes Q = 0.0, 1.0, 1.0, 1.0
        Weights Stokes U = 0.0, 1.0, 1.0, 1.0
        Weights Stokes V = 1.0, 1.0, 1.0, 1.0
        Instrumental profile = 3.0

    [Atmospheres]

        [[Photosphere 2]]
        Name = ph2
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300,
        Reference frame = vertical

            [[[Ranges]]]
            T      = -3000.0, 3000.0
            vmic   = 0.0, 3.0
            v      = -10.0, 10.0
            Bx     = -1000.0, 1000.0
            By     = -1000.0, 1000.0
            Bz     = -1000.0, 1000.0
            ff     = 0.0, 1.0

            [[[Nodes]]]
            T      = 3, 3, 5, 5
            vmic   = 1, 1, 1, 1
            v      = 1, 1, 1, 1
            Bx     = 1, 1, 1, 1
            By     = 1, 1, 1, 1
            Bz     = 1, 1, 1, 1
            ff     = 0, 0, 0, 0

            [[Regularization]]
            T      = None
            vmic   = None
            v      = None
            Bx     = None
            By     = None
            Bz     = None

        [[Chromosphere 1]]
        Name = ch1                                              # Name of the atmosphere component
        Spectral region = spec1                                 # Spectral region to be used for synthesis
        Height = 3.0                                            # Height of the slab
        Line = 10830                                            # 10830, 5876
        Wavelength = 10826, 10833                         # Wavelength range used for synthesis
        Reference atmospheric model = 'chromospheres/model_chromosphere.1d'    #File with model parameters
        Coordinates for magnetic field vector = 'cartesian'
        

            [[[Ranges]]]
            Bx     = -500, 500
            By     = -500, 500
            Bz     = -500, 500
            tau    = 0.1, 2.0
            v      = -10.0, 10.0
            deltav = 3.0, 12.0
            beta   = 1.0, 2.0
            a      = 0.0, 1.0
            ff     = 0.0, 1.0
            

            [[[Nodes]]]
            Bx     = 0, 0, 1, 1
            By     = 0, 0, 1, 1
            Bz     = 0, 0, 1, 1
            tau    = 0, 0, 0, 0
            v      = 0, 0, 0, 0
            deltav = 0, 0, 0, 0
            beta   = 0, 0, 0, 0
            a      = 0, 0, 0, 0
            ff     = 0, 0, 0, 0

        [[Parametric 1]]
        Name = te1
        Spectral region = spec1
        Wavelength = 10826, 10833
        Reference atmospheric model = 'telluric/model_telluric.1d'    # File with model parameters
        Type = Gaussian           # Gaussian, Voigt, MoGaussian, MoVoigt 

            [[[Ranges]]]
            Lambda0 = -1.0, 1.0
            Sigma = 0.3, 0.5
            Depth = 0.2, 0.8
            a = 0.0, 0.2
            ff = 0.0, 1.0
        
            [[[Nodes]]]
            Lambda0 = 0, 0, 0, 0
            Sigma = 0, 0, 0, 0
            Depth = 0, 0, 0, 0
            a = 0, 0, 0, 0
            ff = 0, 0, 0, 0

        [[Straylight 1]]
        Name = st1
        Spectral region = spec1
        Wavelength = 10826, 10833    
        Reference atmospheric model = 'straylight/model_stray.1d'    # File with model parameters

            [[[Ranges]]]
            v = -1.0, 1.0        
            ff = 0.0, 1.0
        
            [[[Nodes]]]
            v = 0, 0, 0, 0        
            ff = 0, 0, 0, 0