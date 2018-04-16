Inversion
===========

Because inversion is often an ill-defined operation, getting a good inversion
contains sometimes a large fraction of inspiration. In order to get a 
satisfactory fit, one usually needs to tune the
number of nodes (for photospheric atmospheres), the weights of the Stokes
parameters, the initial configuration, etc. So don't expect to get 
a good inversion on your first run of Hazel. But if you keep insisting, 
you will get something useful out of the code.

A typical configuration file reads something like in the following. Remember that you have
to adapt it to the specific details of your spectrum.

::

    # Hazel configuration File

    [Working mode]
    Output file = output.h5
    Number of cycles = 1


    [Spectral regions]
        [[Region 1]]
        Name = spec1
        Wavelength = 10826, 10833, 150
        Topology = ph1 -> ch1        
        LOS = 0.0, 0.0, 90.0
        Boundary condition = 1.0, 0.0, 0.0, 0.0       # I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)
        Wavelength file = 'observations/10830.wavelength'
        Wavelength weight file = 'observations/10830.weights'
        Observations file = 'observations/10830_stokes.h5'
        Straylight file = None
        Mask file = None
        Weights Stokes I = 1.0, 1.0, 1.0, 1.0
        Weights Stokes Q = 0.0, 1.0, 1.0, 1.0
        Weights Stokes U = 0.0, 1.0, 1.0, 1.0
        Weights Stokes V = 1.0, 1.0, 1.0, 1.0

    [Atmospheres]

        [[Photosphere 1]]
        Name = ph1
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


This configuration file defines a spectral region in the near infrared that contains a photosphere and a chromosphere, one
after the other. The physical conditions of the photosphere are read from the ``photospheres/model_photosphere.1d``
file. Note that all files in this example are 1D so that the input atmosphere is then shared among
all observed pixels. You might as well use 3D files (as described in :ref:`_input`) with the same number of pixels
as the input, and then a different initial atmosphere can be used for every pixel.

In inversion mode you need to define the observations and the weighting scheme, which is done 
through the ``Wavelength file``, ``Wavelength weight file`` and ``Observations file``, all of them explained
in :ref:`input`. If the files with the observations are 1D, you can simply carry out the inversion
by invoking, assuming that the configuration file is saved on ``conf.ini``:

::

    mod = hazel.Model('conf.ini')
    mod.read_observation()
    mod.open_output()
    mod.invert()
    mod.close_output()

The output is described in :ref:`_output`. For 3D cases, we recommend the user to
use iterators. Otherwise, one needs to manually read the observations at every pixel.
For instance, one can do the inversion for many pixels using 3D atmospheres in serial
model by invoking:

::

    iterator = hazel.iterator(use_mpi=False)    
    mod = hazel.Model('conf.ini')
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

or the following piece of code if many cores are to be used. Remember that
this piece of code needs to be called with a call to the MPI ``mpiexec``

::

    mpiexec -n n_nodes python inversion.py

where ``inversion.py`` contains the following piece of code:

::

    iterator = hazel.iterator(use_mpi=True)
    rank = iterator.get_rank()

    if (rank == 0):    
        mod = hazel.Model('conf_mpi_invh5.ini')
        iterator.use_model(model=mod)
    else:
        iterator.use_model()

    iterator.run_all_pixels()