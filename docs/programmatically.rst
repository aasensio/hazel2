.. _programmatically:
.. include:: hazel_name

Programmatically
================

|hazel2| can be used programmatically, which can be used for simple calculations or
if more advanced synthesis/inversion strategies are required. 

Creating a model
----------------
The first step is to create a model. By default, the model is created in synthesis mode.
The verbosity level can be selected with the ``verbose`` keyword.

::

    import hazel
    mod = hazel.Model(working_mode='synthesis', verbose=True)

The ``Model`` class defines the following two properties: ``atmospheres`` and ``spectrum``. The
first one contains all the atmospheric information, and one interacts with it when changing
and/or accessing the properties of each atmosphere. The second one contains all the spectral
information, including the synthesized Stokes parameters or the observations.

Adding a spectral region
------------------------

The first step is to define a spectral region of interest.
The model class has an ``add_spectral`` method that generates new spectral
regions. A dictionary with the information needs to be passed to the
function. Note that many of them are optional and many are only used
in inversion mode. 

::

    dict = {'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1', 
        'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]}
    mod.add_spectral(dict)
    
The explanation for all keys of the dictionary are the same as those explained in :ref:`configuration`:

* ``Name``
* ``Wavelength`` (optional if a wavelength file is passed)
* ``Topology``
* ``Stokes weights`` (optional)
* ``LOS`` (mandatory for synthesis)
* ``Boundary condition`` (mandatory for synthesis)
* ``Wavelength file`` (optional if a wavelength axis is passed
* ``Wavelength weight file`` (optional)
* ``Observations file`` (optional)
* ``Straylight file`` (optional)
* ``Mask file`` (optional)


Working with chromospheres
--------------------------

Adding a chromosphere
^^^^^^^^^^^^^^^^^^^^^

A new Hazel chromosphere is added by calling the following method of ``Model``:

::

    mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 
        'Line': '10830', 'Wavelength': [10826, 10833]})

* ``Name``
* ``Spectral region``
* ``Height``
* ``Line``
* ``Wavelength``
* ``Reference atmospheric model`` (optional)
* ``Ranges`` (optional)
* ``Nodes`` (optional)

In this case, the properties of the newly created atmosphere can be accessed in the
dictionary element ``mod.atmospheres['ch1']``. All atmospheres are indexed by
the given name.

Changing the parameters of the chromosphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once defined, the model parameters can be changed by calling the ``set_parameters``
method of the atmosphere:

::

    mod.atmospheres['ch1'].set_parameters([0.0,0.0,100.0,1.0,0.0,8.0,1.0,0.0],1.0)

where the parameters are, in order: Bx, By, Bz, :math:`\tau`, v, :math:`\Delta v`, :math:`\beta`, a. The filling
factor is passed as the second parameter.
            
Working with photospheres
--------------------------

Adding a photosphere
^^^^^^^^^^^^^^^^^^^^

A new SIR atmosphere can be added by calling the following method of ``Model``:

::

    mod.add_photosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Spectral lines': [300,], 
        'Wavelength': [10826, 10833]})

* ``Name``
* ``Spectral region``
* ``Spectral lines``
* ``Wavelength``
* ``Reference atmospheric model`` (optional)
* ``Ranges`` (optional)
* ``Nodes`` (optional)

Changing the parameters of the photosphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once defined, the model parameters can be changed by calling the ``set_parameters``
method of the atmosphere:

::

    mod.atmospheres['ch1'].set_parameters(model, ff)

where ``model`` is an array of size (nz,8) with the model atmosphere and ff 
is the filling factor. The columns are, in order: :math:`\log \tau`, T, Pe, vmic, v, Bx, By, Bz.
If Pe is negative, the electron pressure will be calculated in synthesis mode
assuming hydrostatic equilibrium. On the other hand, it is always calculated
in hydrostatic equilibrium in inversion mode.


Working with parametric models
--------------------------

Adding a parametric model
^^^^^^^^^^^^^^^^^^^^

::

    mod.add_parametric({'Name': 'te1', 'Spectral region': 'spec1', 'Type': 'Voigt', 
        'Wavelength': [10826, 10833]})

* ``Name``
* ``Spectral region``
* ``Type``
* ``Wavelength``
* ``Reference atmospheric model`` (optional)
* ``Ranges`` (optional)
* ``Nodes`` (optional)

Changing the parameters of the parametric model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once defined, the model parameters can be changed by calling the ``set_parameters``
method of the atmosphere:

::

    mod.atmospheres['te1'].set_parameters([10833, 0.1, 0.5, 0.0, 1.0])

where the parameters are, in order: :math:`\lambda_0`, :math:`\sigma`, d, a, ff.

Working with straylight
--------------------

Adding a straylight
^^^^^^^^^^^^^^^^^^^^

::

    mod.add_parametric({'Name': 'st1', 'Spectral region': 'spec1',  
        'Wavelength': [10826, 10833]})

* ``Name``
* ``Spectral region``
* ``Wavelength``
* ``Reference atmospheric model`` (optional)
* ``Ranges`` (optional)
* ``Nodes`` (optional)

Changing the parameters of the straylight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once defined, the model parameters can be changed by calling the ``set_parameters``
method of the atmosphere:

::

    mod.atmospheres['st1'].set_parameters([0.0, 1.0])

where the parameters are, in order: v, ff.

Finalizing setup
----------------

Once all spectral regions and atmospheres are added, we need to finalize the
setup by invoking:

::

    mod.setup()


Synthesis
---------

The model can be synthesized by calling the ``synthesize`` method of ``Model``

::

    mod.synthesize()

Finally, the emergent Stokes parameters can be accessed, for each spectral region,
by examining ``mod.spectrum['spec1'].stokes``.
