.. _output:

Output files
===========

Hazel v2.0 can accept several formats for input/output files. 1D formats are not allowed
for output and the output is either HDF5 or FITS. 3D files also save some metadata of 
interest. In the case of HDF5, you can access this metadata using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    version = f.attrs['version']
    date = f.attrs['date']

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

Each specific spectral region defines an HDF5 ``Group``, with two different ``Datasets``. The first
one is ``stokes``, which contains the emergent Stokes profiles for all the cycles and is of size ``(n_pixel,n_cycles,4,n_lambda)``. The
second is the value of the :math:`\chi^2` merit function for each cycle and is of size ``(n_pixel,n_cycles)``.

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
inside the group. Using ``h5py``, they can be easily accessed using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    T = f['ph1']['T']

The shape of the output is ``(n_pixel,n_cycles,nz)`` for each specific parameter of a photospheric model
and ``(n_pixel,n_cycles,1)`` for each parameter of the remaining model atmospheres. For the sake of
clarity, the units of the output are saved as attributes. You can watch the units by invoking:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    print(f['ph1']['T'].attrs['unit'])

FITS files
^^^^^^^^^^
TBD