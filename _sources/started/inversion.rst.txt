.. include:: ../hazel_name

Inversion
===========

Because inversion is often an ill-defined operation, getting a good inversion
contains sometimes a large fraction of inspiration. In order to get a 
satisfactory fit, one usually needs to tune the
number of nodes (for photospheric atmospheres), the weights of the Stokes
parameters, the initial configuration, etc. So don't expect to get 
a good inversion on your first run of Hazel. But if you keep insisting, 
you will get something useful out of the code.

Configuration file
------------------

A typical configuration file reads something like in the following. Remember that you have
to adapt it to the specific details of your spectrum.

::

    # Hazel configuration File

    [Working mode]
    Output file = output.h5
    Number of cycles = 1
    Save all cycles = False


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
all observed pixels. You might as well use 3D files (as described in :ref:`input`) with the same number of pixels
as the input, and then a different initial atmosphere can be used for every pixel.

In inversion mode you need to define the observations and the weighting scheme, which is done 
through the ``Wavelength file``, ``Wavelength weight file`` and ``Observations file``, all of them explained
in :ref:`input`. 

Single inversions
-----------------

Without randomization
^^^^^^^^^^^^^^^^^^^^^

If the files with the observations are 1D, you can simply carry out the inversion
by invoking, assuming that the configuration file is saved on ``conf.ini``:

::

    mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion')
    mod.read_observation()
    mod.open_output()
    mod.invert()
    mod.write_output(
    mod.close_output()

We refer to :ref:`output` for a description of the output. 

With randomization
^^^^^^^^^^^^^^^^^^

Carrying out robust inversions is a difficult task. Randomization is a technique
appropriate to capture the sensitivity of the result to the initial conditions.
In the ideal case, one should obtain the same solution irrespectively of the initial
values of the parameters. However, this will never be the case. The best one can
do is to get the sensitivity of the results to the initial conditions. To this end,
|hazel2| provides a method to randomize the initial conditions. Just pass the maximum
number of randomizations when instantiating the ``Model`` and then you can
carry out inversions by looping over randomizations and saving each result to the
output file. We prefer the user to deal with randomization because of the added
flexibility.

::

    mod = hazel.Model('configurations/conf_single.ini', working_mode='inversion', randomization=2)
    mod.read_observation()
    mod.open_output()
    for loop in range(2):
        mod.invert(randomize=True)
        mod.write_output(randomization=loop)
    mod.close_output()


Multiple inversions
-------------------

For 3D cases, we recommend the user to
use iterators. Otherwise, one needs to manually read the observations at every pixel.

Serial inversions
^^^^^^^^^^^^^^^^^

In serial mode, one can do the inversion of many pixels using 3D atmospheres by using the
following strategy.

Without randomization
"""""""""""""""""""""

::

    iterator = hazel.Iterator(use_mpi=False)    
    mod = hazel.Model('configurations/conf_nonmpi_inv1d.ini', working_mode='inversion', verbose=2)
    iterator.use_model(model=mod)
    iterator.run_all_pixels(start=0)

The first thing to do is to instantiate an ``Iterator`` which will be used to
iterate over all pixels in the file with the observations. The ``Iterator`` class
admits a single keyword ``use_mpi``, which is set to ``False`` by default.
If ``True``, it will use the `Message Passing Interface
<https://en.wikipedia.org/wiki/Message_Passing_Interface/>`_ machinery, for which
you need to have the ``mpi4py`` package installed. If ``False``, it will serially
iterate over all pixels. Then we instantiate the ``Model`` as usual and tell the
``Iterator`` to use this model. Finally, we run the calculation by calling the 
``run_all_pixels`` method of the ``Iterator``. ``run_all_pixels`` admit two parameters:
``start``, set by default to ``0`` that marks the starting point of the inversions and
``end``, set by default to ``None`` (equivalent to using the number of pixels in the file 
with the observations), that mark the last pixel to invert.

With randomization
""""""""""""""""""

The serial mode can also be used with randomization. Because of the complexity, the
iterator takes fully care of randomizations and will simplify your life. It can be
used simply by passing the ``randomization`` keyword during the ``Model`` instantiation.

::

    iterator = hazel.Iterator(use_mpi=False)    
    mod = hazel.Model('configurations/conf_nonmpi_inv1d.ini', working_mode='inversion', verbose=2, randomization=2)
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

Parallel inversions
^^^^^^^^^^^^^^^^^^^

Carrying out inversions of maps in serial mode is probably not the way to go. The
desired strategy is to use computers with several nodes, which can work in parallel.
The specific paralellization in |hazel2| is such that practically a linear scaling
with the number of nodes is sure. It uses a parent-worker architecture, in which a
parent node reads the input and broadcasts it to the nodes, which are then
in charge of carrying out the synthesis/inversions. When finished, the results are
sent back to the parent, who saves the results and sends a new observation to the
available worker.

Without randomization
"""""""""""""""""""""

Running an MPI work needs to be done by using ``mpiexec``, so that you need
to define the following lines in a script:

::

    iterator = hazel.Iterator(use_mpi=True)    
    mod = hazel.Model('configurations/conf_mpi_invh5.ini', working_mode='inversion', rank=iterator.get_rank())
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

The only difference with respect to the serial mode is that you need to pass the
rank of each worker during the instantiation of the ``Model`` because the parent and the
worker will need to do different things with the model. Then just run it with ``mpiexec``, 
defining the number of nodes used in this work:

::

    mpiexec -n n_nodes python inversion.py

Sometimes you need to let OpenMPI or MPICH use all available cores (this is important
if your CPUs are using hyperthreading). To this end, write a file (for example `mf`) with the following
content

::

    localhost max_slots=16

changing the number to the appropriate number of available cores you have. Then run the
code with:


::

    mpiexec -machinefile mf -n n_nodes python inversion.py


With randomization
""""""""""""""""""

Randomization can also be used in parallel. For that, just follow the same approach as before and
pass the ``randomization`` keyword to the ``Model``:

::

    iterator = hazel.Iterator(use_mpi=True)    
    mod = hazel.Model('configurations/conf_mpi_invh5.ini', working_mode='inversion', rank=iterator.get_rank(), randomization=2)
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

