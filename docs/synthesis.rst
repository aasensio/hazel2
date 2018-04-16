Synthesis
===========

This is the simplest use of the code because it should always work once the input 
and configuration files are correctly defined. A typical configuration file for
the synthesis mode follows:

::

    # Hazel configuration File

    [Working mode]
    Output file = output.h5
    

    [Spectral regions]
        [[Region 1]]
        Name = spec1
        Wavelength = 10826, 10833, 150
        Topology = ph1 -> ch1 -> te1 -> st1    
        Stokes weights = 1.0, 1.0, 1.0, 1.0
        LOS = 0.0, 0.0, 90.0
        Boundary condition = 1.0, 0.0, 0.0, 0.0       # I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)
        Wavelength file = 'observations/10830.wavelength'
        Wavelength weight file = 'observations/10830.weights'
    
    [Atmospheres]

        [[Photosphere 1]]
        Name = ph1
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

This configuration file defines a spectral region in the near infrared that contains four different atmospheres, one
after the other. The first one is a photosphere, with the physical conditions been read from the ``photospheres/model_photosphere.1d``
file. Note that all files in this example are 1D so the result is just a single synthesis. But you might
as well use 3D files (as described in :ref:`input`) and the output would be the synthesis for many pixels. The second
component is a chromosphere, the third is used to add a telluric component and the fourth one is
a straylight contamination.

In the case of a 1D model, you can easily do the synthesis by invoking, assuming that the configuration file
is saved on ``conf.ini``:

::

    mod = hazel.Model('conf.ini')
    mod.open_output()
    mod.synthesize()
    mod.close_output()

The output is described in :ref:`output`. For 3D cases, we recommend the user to
use iterators. Otherwise, one needs to manually read the observations at every pixel.
For instance, one can do the synthesis for many pixels using 3D atmospheres in serial
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