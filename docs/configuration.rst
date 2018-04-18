.. _configuration:
.. include:: hazel_name

Configuration
=============

There is a single human-readable configuration file that controls the behavior of the code. It can be
used both for synthesis and inversion.


Example for synthesis
-------

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
    
    [Atmospheres]

        [[Photosphere 2]]
        Name = ph2
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300,
    
        [[Chromosphere 1]]
        Name = ch1                                              # Name of the atmosphere component
        Spectral region = spec1                                 # Spectral region to be used for synthesis
        Height = 3.0                                            # Height of the slab
        Line = 10830                                            # 10830, 5876
        Wavelength = 10826, 10833                         # Wavelength range used for synthesis
        Reference atmospheric model = 'chromospheres/model_chromosphere.1d'    # File with model parameters

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

Example for inversion
-------

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

    [Atmospheres]

        [[Photosphere 2]]
        Name = ph2
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300,

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
        Reference atmospheric model = 'chromospheres/model_chromosphere.1d'    # File with model parameters

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

Working mode
------------

The first part of the configuration file represents very general properties.

::

    [Working mode]
    Output file = output.h5
    Number of cycles = 1
    Maximum iterations = 10
    Relative error = 1e-4

* ``Ouput file``: defines the output file, which is usually an HDF5 or FITS file. It should always be present, otherwise you won't get any output.
* ``Number of cycles`` (optional) : is a global variable to select the number of cycles to carry out during inversion. It can be used to neglect the number of cycles that will be described later in the configuration file.
* ``Maximum iterations`` (optional, default is 10) : maximum number of iterations per cycle to carry out
* ``Relative error`` (optional, default is 1e-4) : relative error when to stop iterating

Spectral regions
----------------

The `spectral regions` are considered as the main objects in |hazel2|. You can
add any number of spectral regions, which will cover the observed region or the
region you desire for your synthesis. One also needs to define which topology of
atmospheres (described below) will produce the synthetic profiles for this region.

::

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

* ``Name``: defines the name of the spectral region. Programmatically, one can have access to the spectrum of each spectral region by using ``mod.spectrum['spec1'].stokes``. This name is also used in the output file to refer to each region.
* ``Wavelength`` (optional) : defines the lower, upper and number of points in the wavelength axis. It can be absent if a file with the wavelength axis is provided.
* ``Topology`` defines the combination of atmospheres that are used to synthesize the Stokes parameters in this spectral region. See :ref:`topology` for more details on the syntax.
* ``LOS`` (mandatory for synthesis) defines the line-of-sight angles: :math:`\theta_\mathrm{LOS}`, :math:`\phi_\mathrm{LOS}` and :math:`\gamma_\mathrm{LOS}`
* ``Boundary condition`` (mandatory for synthesis) defines the boundary condition normalized to the continuum intensity on the quiet Sun at disk center
* ``Wavelength file`` (optional) defines which wavelength file to be used. See :ref:`input` for more information about the format.
* ``Wavelength weight file`` (optional) defines the wavelength weights to be used during inversion, in case one wants to weight parts of the spectrum during the inversion.
* ``Observations file`` (optional) defines the file with the observations. See :ref:`input` for more information.
* ``Straylight file`` (optional) defines the file with the straylight. See :ref:`input` for more information.
* ``Mask file`` (optional) defines a mask to invert only a selection of pixels from an input file. See :ref:`input` for more information.
* ``Weights Stokes`` (optional) defines the weights for all Stokes parameters and cycles. If absent, they will be considered to be 1.

Atmospheres
-----------

The last part of the configuration file defines all atmospheres to be used during the
synthesis. Note that an arbitrary number of atmospheres can be defined. If they are
not used because they are not part of any topology, they will be removed from the
calculation. We distinguish four types of atmospheres:

* Photospheres
* Chromospheres
* Parametric
* Straylight

Photospheres
^^^^^^^^^^^^

Photospheres are always in the lower part of the atmosphere and so need to be defined
in the first level of the topology. They are synthesized in local thermodynamic equilibrium
using SIR. The magnetic field is always given with respect to the line-of-sight (contrary
to those of chromospheres).

::

    [Atmospheres]

        [[Photosphere 1]]
        Name = ph1
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300, 301

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


* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` : defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time.
* ``Spectral lines`` : it is a comma-separated list of lines to synthesize from the :ref:`photospheric_lines`. Note that if you only want one line, you should use a comma at the end. The list of available lines
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : not yet implemented

Chromospheres
^^^^^^^^^^^^^

Chromospheres are synthesized with slabs of constant physical properties always above photospheres (if any). An arbitrary
number of chromospheres can be used, either with filling factor or one above the other. Note that the magnetic
field is always given with respect to the local vertical, contrary to photospheres, which are given 
relative to the line of sight.

::

    [Atmospheres]

        [[Chromosphere 1]]
        Name = ch1                                              # Name of the atmosphere component
        Spectral region = spec1                                 # Spectral region to be used for synthesis
        Height = 3.0                                            # Height of the slab
        Line = 10830                                            # 10830, 5876
        Wavelength = 10826, 10833                         # Wavelength range used for synthesis
        Reference atmospheric model = 'chromospheres/model_chromosphere.1d'    # File with model parameters

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

* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` : defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time.
* ``Line`` : which of the He I lines to consider (5876, 10830, ...)
* ``Height`` : height of the slab in arcsec.
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : not yet implemented

Parametric
^^^^^^^^^^

Parametric atmospheres are used to synthesize any systematics that can be affecting the
observations. Things like telluric lines, fringes, smooth continua form part of this. You can
always correct from them during data reduction, but |hazel2| allows you to add them
during the fit. The curently available parametric atmosphere is just a Voigt function.

::

    [Atmospheres]

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

* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` : defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time.
* ``Type`` : type of parametric atmosphere, from the available selection ``Voigt``/``MoVoigt``
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : not yet implemented

Straylight
^^^^^^^^^^

Straylight components are always added to the final spectrum with a filling factor and a possible velocity shift.
::

    [Atmospheres]

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

* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` : defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time.
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : not yet implemented