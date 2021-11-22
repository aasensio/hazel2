.. include:: ../hazel_name
.. _configuration:


Step-by-step description
========================

All configuration options are described here. Note that some of them are optional
and can be absent from the configuration file. If this is the case, a default
value will be used.

Working mode
------------

The first part of the configuration file represents very general properties.

::

    [Working mode]
    Output file = output.h5
    Number of cycles = 1
    Maximum iterations = 10
    Relative error = 1e-4
    Backtracking = 'brent'
    Save all cycles = False

* ``Ouput file``: defines the output file, which is usually an HDF5 or FITS file. It should always be present, otherwise you won't get any output.
* ``Number of cycles`` (optional) : is a global variable to select the number of cycles to carry out during inversion. It can be used to neglect the number of cycles that will be described later in the configuration file.
* ``Maximum iterations`` (optional, default is 10) : maximum number of iterations per cycle to carry out
* ``Relative error`` (optional, default is 1e-4) : relative error when to stop iterating
* ``Backtracking`` (optional, default is ``brent``) : method to be used for the computation of the optimal Levenberg-Marquardt parameter. Two options are available: ``parabolic`` (which uses parabolic interpolation) and ``brent`` (which uses the Brent algorithm, which should be more efficient)
* ``Save all cycles`` (optional, default is ``False``) : True if you want to save the result of all cycles in the inversion.

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
        Instrumental profile = 3.0

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
* ``Instrumental profile`` (optional) defines the instrumental profile. It can be absent, equal to ``None``, a float giving the width of a Gaussian PSF or a file with the PSF (given with two columns with wavelength displacement in A and PSF).

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
using SIR. The magnetic field is *always* given with respect to the line-of-sight, like
in SIR (contrary to those of chromospheres) and the reference direction for Q is
then defined as the zero azimuth (like in SIR).

::

    [Atmospheres]

        [[Photosphere 1]]
        Name = ph1
        Reference atmospheric model = 'photospheres/model_photosphere.1d'
        Spectral region = spec1
        Wavelength = 10826, 10833
        Spectral lines = 300, 
        NLTE = False
        Temperature change to recompute departure coefficients = 5.0
        Reference frame = line-of-sight

            [[[Ranges]]]
            T      = -3000.0, 3000.0
            vmic   = 0.0, 3.0
            v      = -10.0, 10.0
            Bx     = -1000.0, 1000.0
            By     = -1000.0, 1000.0
            Bz     = -1000.0, 1000.0
            ff     = 0.0, 1.0
            vmac   = 0.0, 5.0

            [[[Nodes]]]
            T      = 3, 3, 5, 5
            vmic   = 1, 1, 1, 1
            v      = 1, 1, 1, 1
            Bx     = 1, 1, 1, 1
            By     = 1, 1, 1, 1
            Bz     = 1, 1, 1, 1
            ff     = 0, 0, 0, 0
            vmac   = 0, 0, 0, 0

            [[[Regularization]]]
            T      = None
            vmic   = None
            v      = None
            Bx     = None
            By     = None
            Bz     = None
            vmac   = None


* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` (optional): defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time. If absent or `None', then the whole spectral region is synthesized.
* ``Spectral lines`` : it is a comma-separated list of lines to synthesize from the :ref:`photospheric_lines`. Note that if you only want one line, you should use a comma at the end. The list of available lines
* ``NLTE`` (optional, default is ``False``): defines whether the line should be treated in NLTE. At the moment, only valid for the Ca II 8542 A line (line number 301).
* ``Temperature change to recompute departure coefficients`` (optional, default is 0): temperature change during the inversion iterations to recompute the departure coefficients.
* ``Reference frame``: reference frame to give the components of the magnetic field: ``line-of-sight`` or ``vertical`` (default).
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range using a logit transform (with a small :math:`\epsilon` to avoid under/overflow when close to the border).
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions. It is possible to couple parameters with those of another atmosphere by using the name of the atmosphere instead of the number of nodes.
* ``Regularization`` : add regularization to the parameters. See :ref:`regularization`.

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

* ``Name`` : defines the name of the atmosphere. This will be used in the output file to refer to the parameters of this specific atmosphere.
* ``Reference atmospheric model`` (optional) : defines the input file for this atmosphere. The format is described in :ref:`input`. If the format is 1D, it will be used for all pixels (in inversion mode). If you want a different model for all pixels, use 3D formats.
* ``Spectral region`` : defines the spectral region associated with this atmosphere.
* ``Wavelength`` (optional): defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time. If absent or `None', then the whole spectral region is synthesized.
* ``Line`` : which of the He I lines to consider (5876, 10830, ...)
* ``Reference frame`` : it defines the reference system in which the magnetic field is measured. ``line-of-sight`` or ``vertical`` (traditional mode for |hazel2|). If absent, ``vertical`` is used by default.
* ``Height`` : height of the slab in arcsec.
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions. It is possible to couple parameters with those of another atmosphere by using the name of the atmosphere instead of the number of nodes.
* ``Regularization`` : add regularization to the parameters. See :ref:`regularization`.
* ``Coordinates for magnetic field vector`` (optional) : defines the coordinates used for defining the magnetic field vector. If ``cartesian`` is used, then you need to define ``Bx``, ``By`` and ``Bz`` as variables in the ranges, nodes, etc. If ``spherical`` is used, then the variables are termed ``B``, ``thB`` and ``phiB``.

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
* ``Wavelength`` (optional): defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time. If absent or `None', then the whole spectral region is synthesized.
* ``Type`` : type of parametric atmosphere, from the available selection ``Voigt``/``MoVoigt``
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : add regularization to the parameters. See :ref:`regularization`.

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
* ``Wavelength`` (optional): defines the ranges to be used for the synthesis of this atmosphere. This is interesting if you only want this atmosphere to synthesize part of the observed spectrum, which will affect the computing time. If absent or `None', then the whole spectral region is synthesized.
* ``Ranges`` : ranges of variation of each parameter. If ``None``, consider it unconstrained. If not, it will be constrained to the range.
* ``Nodes`` : defines the number of nodes in each cycle when doing inversions
* ``Regularization`` : add regularization to the parameters. See :ref:`regularization`.
