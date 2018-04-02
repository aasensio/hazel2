.. inputOutput-label:
Input/output files
===================

Both input and output files for Hazel are NetCDF files.

Input files
-----------

The input file constains the observations and information about the
observing position and boundary condition. The file consists of the
following variables:

-  lambda: vector of size *nlambda* containing the wavelength axis with
   respect to the center of the multiplet (10829.0911 Angstrom for the multiplet at 10830 Angstrom).

-  map: array of size *(npixel,8,nlambda)* containing the Stokes vector
   :math:`(I,Q,U,V)` and the associated standard deviation of the noise
   :math:`(\sigma_I,\sigma_Q,\sigma_U,\sigma_V)`.

-  boundary: array of size *(npixel,4)* containing the boundary
   condition for every inverted pixel.

-  height: vector of size *npixel* which contains the height of the
   slabs for every pixel.

-  obs\_theta: vector of size *npixel* which contains the observing
   angle :math:`\theta` for every pixel.

-  obs\_gamma: vector of size *npixel* which contains the observing
   angle :math:`\gamma` that defines the positive reference for Stokes
   :math:`Q` for every pixel.

-  mask: array of size *nx,ny* which tells whether this pixel will be
   inverted.

-  normalization: variable indicating whether the profiles are
   normalized to the peak amplitude or the continuum of Stokes
   :math:`I`.

-  pars: array of size *npixel,npars* which contains the initial value
   for the model parameters. These will be used to reinvert some pixels
   or, for instance, to refine the ambiguous solutions.

The routine ``gen_netcdf.pro`` on the directory ``IDL_routines`` and the
``genNetCDF.py`` on ``pyRoutines`` shows functions that generate such a
file by passing all the variables as parameters. The order of pars is
the following, depending on the number of slabs:

-  1-component (vector of size 8): :math:`B`, :math:`\theta_B`,
   :math:`\chi_B`, :math:`\tau`, :math:`v_\mathrm{dop}`, :math:`a`,
   :math:`v_\mathrm{mac}`, :math:`\beta`

-  2-component 1+1 with same field (vector of size 11): :math:`B`,
   :math:`\theta_B`, :math:`\chi_B`, :math:`\tau_1`, :math:`\tau_2`,
   :math:`v_\mathrm{dop}`, :math:`a`, :math:`v_\mathrm{mac1}`,
   :math:`v_\mathrm{mac2}`, :math:`\beta`, :math:`\beta_2`

-  2-component 1+1 with different field (vector of size 15):
   :math:`B_1`, :math:`\theta_{B1}`, :math:`\chi_{B1}`, :math:`B_2`,
   :math:`\theta_{B2}`, :math:`\chi_{B2}`, :math:`\tau_1`,
   :math:`\tau_2`, :math:`v_\mathrm{dop}`, :math:`v_\mathrm{dop2}`,
   :math:`a`, :math:`v_\mathrm{mac1}`, :math:`v_\mathrm{mac2}`,
   :math:`\beta`, :math:`\beta_2`

-  2-component 2 with different field with filling factor (vector of
   size 16): :math:`B_1`, :math:`\theta_{B1}`, :math:`\chi_{B1}`,
   :math:`B_2`, :math:`\theta_{B2}`, :math:`\chi_{B2}`, :math:`\tau_1`,
   :math:`\tau_2`, :math:`v_\mathrm{dop}`, :math:`v_\mathrm{dop2}`,
   :math:`a`, :math:`v_\mathrm{mac1}`, :math:`v_\mathrm{mac2}`,
   :math:`\mathrm{ff}`, :math:`\beta`, :math:`\beta_2`

Data preparation
----------------
Some data preprocessing has to be done in order to have reliable inversions with
Hazel. The process takes the following steps:

-  **Data normalization**: the data has to be normalized to the local continuum. This
   is sometimes slightly difficult because the nearby Si I line has strong wings and one should
   use that pseudocontinuum. The very first step would be to remove large scale variations
   of the continuum, so that it is as flat as possible (perhaps removing fringes if you have
   any). Then, you proceed to remove the influence
   of the Si I line. What our people typically use is to fit the Si I line using
   almost all its blue wing and only part of the red wing. I guess you can do it using
   an inversion code like SIR or use a Voigt function. You probably want to get photospheric
   information from your observations, so maybe SIR is a better option. Once you have
   the Si I line fitted, just extend the synthetic wing towards the He I line and then normalize
   the spectrum by the synthetic profile. This way, you'll have the He I triplet correctly
   normalized.
   If the data is off-limb, things are typically easier because there is no continuum but
   sometimes there is some stray-light that can give you a headache. In this case, 
   the input should be normalized by the peak emission.

-  **Wavelength calibration**: the data has to be wavelength calibrated. How to do it depends 
   on whether you want an absolute calibration of velocities or not. If you want such absolute scale, the best is to do the wavelength
   calibration using telluric lines and then transform everything to the Sun using the relative
   velocity between the observed region and the Earth. If not, maybe using some weak 
   surrounding photospheric lines is enough. Note that all wavelengths are given with respect
   to the center of the multiplet, which is 10829.0911 Angstrom for the 10830 Angstrom one)

-  **Computation of the boundary condition and heliocentric angle**: every pixel should be labeled with its heliocentric angle (this is important for
   observations close to the limb, where mu is changing fast) and its boundary condition.
   So, you need to get a map of heliocentric angles together with your map of observed Stokes
   profiles. Concerning the boundary condition, it is enough to compute the ratio between
   the continuum intensity at every pixel and the average at the same heliocentric angle.

-  **Rotation of the reference system**: it is important to understand which is the reference direction for positive Stokes
   Q in the observations. Note that the output of the code depends on the :math:`\gamma` angle, which exactly
   defines this positive Q direction. Two possibilities appear. The first one is to set :math:`\gamma` in the code
   so that you understand which is the reference direction in the code and then rotate the Stokes Q and U data
   so that the reference direction for Stokes Q is aligned with that of the code. The second possibility is to
   keep the data as it is and then put the appropriate value of :math:`\gamma` in the code to make both reference
   directions equal. This is usually not difficult, but it requires to understand which is the reference direction
   for the telescope, which is sometimes difficult to get. It is always a good advice to have the scattering
   geometry in mind and try to adapt it to your observations. See :ref:`image_geometry` for more information.


Output files
------------

The results of the inversion are saved on two files defined on the
configuration file. The file with the inverted
profiles contains the following variables:

-  lambda: vector of size *nlambda* containing the wavelength axis with
   respect to the center of the multiplet (10829.0911 Angstrom for the multiplet at 10830 Angstrom)

-  map: array of size *(npixel,4,nlambda)* containing the synthetic
   Stokes vector :math:`(I,Q,U,V)` for every pixel.

The file with the inverted parameters contains the following variable:

-  map: array of size *(npixel,ncolumns)* containing the parameters of
   the inversion.

The number of columns depends on the selected model:

-  One-slab: nine columns with the vector
   :math:`(B,\theta_B,\chi_B,h,\tau,v_\mathrm{th},a,v_\mathrm{mac},\beta)`.

-  Two-slab with same magnetic field: eleven columns with the vector
   :math:`(B,\theta_B,\chi_B,h,[\tau]_1,[\tau]_2,v_\mathrm{th},a,[v_\mathrm{mac}]_1,
   [v_\mathrm{mac}]_2,\beta,\beta_2)`.

-  Two-slab with different magnetic field: fifteen columns with the
   vector
   :math:`([B]_1,[\theta_B]_1,[\chi_B]_1,[B]_2,[\theta_B]_2,[\chi_B]_2,
   h,[\tau]_1,[\tau]_2,[v_\mathrm{th}]_1,[v_\mathrm{th}]_2,a,[v_\mathrm{mac}]_1,[v_\mathrm{mac}]_2,\beta,\beta_2)`.

The file ``read_results.pro`` on the ``RunMPI`` directory shows how to
read the files from IDL.