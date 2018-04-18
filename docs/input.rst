.. _input:
.. include:: hazel_name

Input files
===========

|hazel2| can accept several formats for input/output files. 1D formats are specific
for Hazel, but 3D formats are defined in standard HDF5 or FITS formats.

Wavelength files
----------------

Wavelength
^^^^^^^^^^
Wavelength files are defined with a text file with a single column and a header, and
defines the wavelength axis in Angstrom of each spectral region.

::
    
    # lambda
    1.082600000000000000e+04
    1.082604697986577230e+04
    1.082609395973154278e+04
    1.082614093959731508e+04
    1.082618791946308738e+04
    1.082623489932885968e+04
    ...

Wavelength weights
^^^^^^^^^^^^^^^^^^
Each wavelength point can be weighted differently in the four Stokes parameters. To this end, you can
provide a file with the weights in the following format.

::
    
    # WeightI WeightQ WeightU WeightV
       1.0     1.0      1.0    1.0
       1.0     1.0      1.0    1.0
       1.0     1.0      1.0    1.0
       1.0     1.0      1.0    1.0
       ......


Observations files
------------------

We describe now the files that contain the observations needed for the inversion mode. Observations
can be described in 1D files for single-pixel inversions, or in 3D files with different flavors for
many-pixels inversions.

1D files
^^^^^^^^

These are text files giving the necessary information for the inversion. The first data line
is the line-of-sight angles :math:`\theta_\mathrm{LOS}`, :math:`\phi_\mathrm{LOS}` and :math:`\gamma_\mathrm{LOS}`.
The second line gives the boundary condition in units of the continuum intensity of the quiet Sun at disk center.
Finally, one has to list the value of the Stokes parameters and the standard deviation of the
noise for each wavelength and Stokes parameter. An example follows:

::

    # LOS theta_LOS, phi_LOS, gamma_LOS
    0.0 0.0 90.0

    # Boundary condition I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)
    1.0 0.0 0.0 0.0

    # SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV
    9.398283088919315853e-01 2.307830414199267630e-04 -5.121676330588738812e-05 1.457157835263802345e-04 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04
    9.354598051184511709e-01 1.898170981935559632e-04 -1.157160550018296303e-04 -8.932093956208153021e-05 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04
    9.304626067718699822e-01 3.837834915890673399e-05 1.166320038345326575e-04 1.087068643459882281e-04 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04 1.000000000000000048e-04
    ....


HDF5 3D files
^^^^^^^^^^^^^

HDF5 files with observations are defined with four double-precision datasets: ``stokes``, ``sigma``, ``LOS`` and ``boundary``.
The first two and the boundary conditions have size ``(n_pixel,n_lambda,4)``, 
and contain the four Stokes parameters and standard deviation of the noises for all pixels, together with the
boundary conditions for each pixel. The LOS is of size ``(n_pixel,3)`` and contains the 
:math:`\theta_\mathrm{LOS}`, :math:`\phi_\mathrm{LOS}` and :math:`\gamma_\mathrm{LOS}` angles for
each pixel. In the following we show how to
create a sample file:

::

    n_pixel = 100
    n_lambda = 150

    stokes_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
    sigma_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64))
    los_3d = np.zeros((n_pixel,3), dtype=np.float64)
    boundary_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

    f = h5py.File('observations/10830_stokes.h5', 'w')
    db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
    db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
    db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
    db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
    db_stokes[:] = stokes_3d
    db_sigma[:] = sigma_3d
    db_los[:] = los_3d
    db_boundary[:] = boundary_3d
    f.close()

FITS 3D files
^^^^^^^^^^^^^
TBD.

Mask file
^^^^^^^^^
Sometimes it is interesting to reinvert some pixels of a map that, for some reason, did not
go well in a previous inversion. Instead of rewriting a new file with these pixels, Hazel2
allows you to use a file with a mask. This mask marks with 1 those pixels that need to be
inverted and with zero those that will not be inverted. For the moment, we do not have
the option of merging a previous inversion with a new one, so that the output file
will be empty in those pixels that are not inverted. You will need to merge the
two output files yourself. Since this mask is appropriate only for 3D inversions, 
it needs to be in one of the 3D formats. There is a single dataset ``mask`` of
type ``np.int8``. For the moment, only 0 and 1 are used, but we prefer to use
an integer in case more elaborate masks can be used in the future. If the file is
``None`` or the keyword is absent, all pixels will be inverted.

::

    n_pixel = 100    
    mask_3d = np.zeros((n_pixel,), dtype=np.int8)

    f = h5py.File('observations/10830_mask.h5', 'w')
    db_mask = f.create_dataset('mask', mask_3d.shape, dtype=np.int8)    
    db_mask[:] = mask_3d    
    f.close()

Photospheric models
-------------------

Photospheric models are used both in synthesis and inversion. In synthesis mode, they are used to
obtain the emergent Stokes parameters. In inversion mode, they are used as an initial reference model that will
be perturbed with nodes until a fit is obtained for the observed Stokes parameters.

1D files
^^^^^^^^

These are text files which tabulates the depth dependence as a function of the log optical depth at 500 nm
of the temperature [K], electron pressure [cgs], 
microturbulent velocity [cm/s], bulk velocity [cm/s], and the cartesian components of the magnetic field,
Bx, By and Bz [G]. Additionally, the filling factor is given in the header. An example follows:

::

    ff
    1.0

    logtau     T        Pe           vmic        v            Bx           By         Bz
    1.2000   8879.7  2.99831E+03  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00   
    1.1000   8720.2  2.46927E+03  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00   
    1.0000   8551.0  1.98933E+03  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00   
    0.9000   8372.2  1.56782E+03  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00   
    0.8000   8183.7  1.20874E+03  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00   
    0.7000   7985.6  9.11633E+02  0.000E+00  0.0000E+00   5.0000E+02    0.0000E+00  0.0000E+00
    ...

HDF5 3D files
^^^^^^^^^^^^^

HDF5 files with model photospheres are defined with two double-precision datasets: ``model`` and ``ff``. The first
one has size ``(n_pixel,nz,8)``, containing the depth dependence of the 8 variables for all pixels. The second 
one has size ``(n_pixel,)``, containing the filling factor for each pixel. In the following we show how to
create a sample file:

::

    n_pixel = 100
    nz = 50
    
    model_3d = np.zeros((n_pixel,nz,8), dtype=np.float64)
    ff_3d = np.zeros((n_pixel,), dtype=np.float64)

    f = h5py.File('photospheres/model_photosphere.h5', 'w')
    db_model = f.create_dataset('model', model_3d.shape, dtype=np.float64)
    db_ff = f.create_dataset('ff', ff_3d.shape, dtype=np.float64)
    db_model[:] = model_3d
    db_ff[:] = ff_3d
    f.close()

FITS 3D files
^^^^^^^^^^^^^
TBD

Chromospheric models
--------------------

Chromospheric models are used both in synthesis and inversion. In synthesis mode, they are used to
obtain the emergent Stokes parameters. In inversion mode, they are used as an initial reference model that will
be perturbed with nodes until a fit is obtained for the observed Stokes parameters.

1D files
^^^^^^^^

These are text files which tabulates the three cartesian components of the magnetic field [G], 
the optical depth of the slab, the bulk velocity [km/s], the Doppler width of the line [km/s], 
the enhancement factor beta, the damping a and the filling factor. An example follows:

::

    Bx [G]   By [G]   Bz [G]   tau    v [km/s]     deltav [km/s]   beta    a     ff
    0.0       0.0     0.0      1.0      0.0             8.0        1.0    0.0    1.0


HDF5 3D files
^^^^^^^^^^^^^

HDF5 files with model chromospheres are defined with two double-precision datasets: ``model`` and ``ff``. The first
one has size ``(n_pixel,8)``, containing the depth dependence of the 8 variables for all pixels. The second 
one has size ``(n_pixel,)``, containing the filling factor for each pixel. In the following we show how to
create a sample file:

::

    n_pixel = 100
    
    model_3d = np.zeros((n_pixel,8), dtype=np.float64)
    ff_3d = np.zeros((n_pixel,), dtype=np.float64)

    f = h5py.File('photospheres/model_chromosphere.h5', 'w')
    db_model = f.create_dataset('model', model_3d.shape, dtype=np.float64)
    db_ff = f.create_dataset('ff', ff_3d.shape, dtype=np.float64)
    db_model[:] = model_3d
    db_ff[:] = ff_3d
    f.close()

Parametric models
-----------------
TBC

Straylight models
-------------------

Straylight models contain a velocity [in km/s] and a filling factor. Additionally, one needs to define the 
intensity spectrum of the stray light contamination, which is given in units of the continuum intensity at
disk center.

1D files
^^^^^^^^

::

    v   ff
    0.0  1.0

    Stokes I
    0.0
    0.0
    0.0
    1.0
    1.0
    1.0
    1.0
    ....


HDF5 3D files
^^^^^^^^^^^^^

HDF5 files with model straylight contaminations are defined with three double-precision datasets: ``profile``, ``model`` and ``ff``. The first
one has size ``(n_pixel,n_wavelength)`` (which needs to conform with the number of wavelength of the spectral region
associated to this atmosphere) and contains the intensity spectrum in units of the continuum intensity
at disk center. The second one has size ``(n_pixel,8)`` and contains the model variable, which in this case
reduces just to the velocity shift. The last one
has size ``(n_pixel,)`` and contains the filling factor for each pixel. In the following we show how to
create a sample file:

::

    n_pixel = 100
    
    model_3d = np.zeros((n_pixel,8), dtype=np.float64)
    ff_3d = np.zeros((n_pixel,), dtype=np.float64)
    profile_3d = np.zeros((n_pixel,150), dtype=np.float64)

    f = h5py.File('straylight/model_straylight.h5', 'w')
    db_model = f.create_dataset('model', model_3d.shape, dtype=np.float64)
    db_profile = f.create_dataset('profile', profile_3d.shape, dtype=np.float64)
    db_ff = f.create_dataset('ff', ff_3d.shape, dtype=np.float64)
    db_model[:] = model_3d
    db_profile[:] = profile_3d
    db_ff[:] = ff_3d
    f.close()