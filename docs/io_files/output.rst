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
    T_err = f['ph1']['T_err']
    T_nodes = f['ph1']['T_nodes']

The shape of the output is ``(n_pixel,n_randomizations,n_cycles,nz)`` for each specific parameter of a photospheric model
and ``(n_pixel,n_randomizations,n_cycles,1)`` for each parameter of the remaining model atmospheres. The location of the nodes
for each variable has dimensions ``(n_pixel,n_randomizations,n_cycles)'' and return an array with the value of log :math:`\tau`
for each node. Additionally, the error associated with each node can also be found. For the sake of
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

Finally, you can draw the tree of the output file, which surely clarifies the format by
invoking ``hazel.util.show_tree('output.h5')``:

::

    output.h5
        ├─ ch1
        │  ├─ B -> (1, 1, 1, 1)  float64
        │  ├─ B_err -> (1, 1, 1)  object
        │  ├─ B_nodes -> (1, 1, 1)  object
        │  ├─ Bx -> (1, 1, 1, 1)  float64
        │  ├─ Bx_err -> (1, 1, 1)  object
        │  ├─ Bx_nodes -> (1, 1, 1)  object
        │  ├─ By -> (1, 1, 1, 1)  float64
        │  ├─ By_err -> (1, 1, 1)  object
        │  ├─ By_nodes -> (1, 1, 1)  object
        │  ├─ Bz -> (1, 1, 1, 1)  float64
        │  ├─ Bz_err -> (1, 1, 1)  object
        │  ├─ Bz_nodes -> (1, 1, 1)  object
        │  ├─ a -> (1, 1, 1, 1)  float64
        │  ├─ a_err -> (1, 1, 1)  object
        │  ├─ a_nodes -> (1, 1, 1)  object
        │  ├─ beta -> (1, 1, 1, 1)  float64
        │  ├─ beta_err -> (1, 1, 1)  object
        │  ├─ beta_nodes -> (1, 1, 1)  object
        │  ├─ deltav -> (1, 1, 1, 1)  float64
        │  ├─ deltav_err -> (1, 1, 1)  object
        │  ├─ deltav_nodes -> (1, 1, 1)  object
        │  ├─ ff -> (1, 1, 1, 1)  float64
        │  ├─ ff_err -> (1, 1, 1)  object
        │  ├─ ff_nodes -> (1, 1, 1)  object
        │  ├─ phiB -> (1, 1, 1, 1)  float64
        │  ├─ phiB_err -> (1, 1, 1)  object
        │  ├─ phiB_nodes -> (1, 1, 1)  object
        │  ├─ tau -> (1, 1, 1, 1)  float64
        │  ├─ tau_err -> (1, 1, 1)  object
        │  ├─ tau_nodes -> (1, 1, 1)  object
        │  ├─ thB -> (1, 1, 1, 1)  float64
        │  ├─ thB_err -> (1, 1, 1)  object
        │  ├─ thB_nodes -> (1, 1, 1)  object
        │  ├─ v -> (1, 1, 1, 1)  float64
        │  ├─ v_err -> (1, 1, 1)  object
        │  └─ v_nodes -> (1, 1, 1)  object
        ├─ ph2
        │  ├─ Bx -> (1, 1, 1, 73)  float64
        │  ├─ Bx_err -> (1, 1, 1)  object
        │  ├─ Bx_nodes -> (1, 1, 1)  object
        │  ├─ By -> (1, 1, 1, 73)  float64
        │  ├─ By_err -> (1, 1, 1)  object
        │  ├─ By_nodes -> (1, 1, 1)  object
        │  ├─ Bz -> (1, 1, 1, 73)  float64
        │  ├─ Bz_err -> (1, 1, 1)  object
        │  ├─ Bz_nodes -> (1, 1, 1)  object
        │  ├─ T -> (1, 1, 1, 73)  float64
        │  ├─ T_err -> (1, 1, 1)  object
        │  ├─ T_nodes -> (1, 1, 1)  object
        │  ├─ ff -> (1, 1, 1, 1)  float64
        │  ├─ ff_err -> (1, 1, 1)  object
        │  ├─ ff_nodes -> (1, 1, 1)  object
        │  ├─ vmac -> (1, 1, 1, 1)  float64
        │  ├─ vmac_err -> (1, 1, 1)  object
        │  ├─ vmac_nodes -> (1, 1, 1)  object
        │  ├─ log_tau -> (73,)  float64
        │  ├─ v -> (1, 1, 1, 73)  float64
        │  ├─ v_err -> (1, 1, 1)  object
        │  ├─ v_nodes -> (1, 1, 1)  object
        │  ├─ vmic -> (1, 1, 1, 73)  float64
        │  ├─ vmic_err -> (1, 1, 1)  object
        │  └─ vmic_nodes -> (1, 1, 1)  object
        └─ spec1
            ├─ aic -> (1, 1, 1)  float64
            ├─ bic -> (1, 1, 1)  float64
            ├─ chi2 -> (1, 1, 1)  float64
            ├─ stokes -> (1, 1, 1, 4, 164)  float64
            └─ wavelength -> (164,)  float64


FITS files
^^^^^^^^^^
TBD