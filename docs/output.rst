.. _output:

Output files
===========

Hazel v2.0 can accept several formats for input/output files. 1D formats are not allowed
for output and the output is either HDF5 or FITS.

Synthesis mode
--------------

HDF5 files
^^^^^^^^^^

In synthesis mode, all spectral regions are saved as HDF5 ``Datasets``. Using ``h5py``, they
can be easily accessed using:

::

    import h5py
    f = h5py.File('output.h5', 'r')
    f['spec1'].shape

The shape of the output is (n_pixel,4,n_lambda) for each specific spectral region.

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

The shape of the output is (n_pixel,n_cycles,nz) for each specific parameter of a photospheric model
and (n_pixel,n_cycles,1) for each parameter of the remaining model atmospheres.

FITS files
^^^^^^^^^^
TBD