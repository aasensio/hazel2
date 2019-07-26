.. include:: ../hazel_name
.. _output:


Output files
============

|hazel2| can accept several formats for input/output files. 1D formats are not allowed
for output and the output is either HDF5 or FITS. 3D files also save some metadata of 
interest. In the case of HDF5, you can access this metadata using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    version = f.attrs['version']
    date = f.attrs['date']

HDF5 files contain some attributes that can help the user reconstruct the
specific options used for a certain run of |hazel2|. The following are of interest:

::

    print(f.attrs['version'])
    print(f.attrs['date'])
    print(f.attrs['configuration'])

The first one defines the version of |hazel2| used. The second one defines the date of creation of the file.
The last one is a copy of the configuration file used for running the code.

Synthesis mode
--------------

HDF5 files
^^^^^^^^^^

In synthesis mode, all spectral regions are saved as HDF5 ``Datasets``. Using ``h5py``, they
can be easily accessed using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    stokes = f['spec1']['stokes']
    chi2 = f['spec1']['chi2']
    aic = f['spec1']['aic']
    bic = f['spec1']['bic']

Each specific spectral region defines an HDF5 ``Group``, with two different ``Datasets``. The first
one is ``stokes``, which contains the emergent Stokes profiles for all the cycles and is of size ``(n_pixel,n_randomizations,n_cycles,4,n_lambda)``. The
second is the value of the :math:`\chi^2` merit function for each cycle and is of size ``(n_pixel,n_cycles)``.
Additionally, we output the value of the `Akaike Information Criterion (AIC) <https://en.wikipedia.org/wiki/Akaike_information_criterion>`_
and the `Bayesian Information Criterion (BIC) <https://en.wikipedia.org/wiki/Bayesian_information_criterion>`_ that can
be useful for model comparison.
If you do not remember the ordering of the indices, you can check them from the file directly by invoking:

::

    for i in range(len(f['spec1']['stokes'].dims):
        print(f['spec1']['stokes'].dims[i].label)


FITS files
^^^^^^^^^^
TBD

Inversion mode
--------------

HDF5 files
^^^^^^^^^^

In inversion mode, apart from the synthetic Stokes profiles obtained for each spectral region,
the model parameters of the active atmospheres are also included in the output. Each atmosphere
defines an HDF5 ``Group``, and each variable of each atmosphere is then an HDF5 ``Datasets``
inside the group. The group defines some attributes that allow the user to understand
the inversion procedure followed by the code. Using ``h5py``, they can be easily accessed using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    T = f['ph1']['T']

The shape of the output is ``(n_pixel,n_randomizations,n_cycles,nz)`` for each specific parameter of a photospheric model
and ``(n_pixel,n_randomizations,n_cycles,1)`` for each parameter of the remaining model atmospheres. For the sake of
clarity, the units of the output are saved as attributes. You can watch the units by invoking:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    print(f['ph1']['T'].attrs['unit'])

If you do not remember the ordering of the indices, you can check them from the file directly by invoking:

::

    for i in range(len(f['ph1']['T'].dims):
        print(f['ph1']['T'].dims[i].label)

The number of cycles, the ranges and the regularization parameters can be obtained as a dictionary using:

::

    import ast
    cycles = ast.literal_eval(f['ph1'].attrs['cycles'])
    regularization = ast.literal_eval(f['ph1'].attrs['regularization'])
    ranges = ast.literal_eval(f['ph1'].attrs['ranges'])

FITS files
^^^^^^^^^^
TBD