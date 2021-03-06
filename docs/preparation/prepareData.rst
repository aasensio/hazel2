.. include:: ../hazel_name

Data preparation		
-----------------		

Some data preprocessing has to be done in order to have reliable inversions with		
Hazel. The process takes the following steps:		
		
-  **Data normalization**: the data has to be normalized to the local quiet Sun continuum. This		
   is sometimes slightly difficult because the nearby Si I line has strong wings and one should		
   use that pseudocontinuum. The very first step would be to remove large scale variations		
   of the continuum, so that it is as flat as possible (perhaps removing fringes if you have		
   any). Then, if you are only inverting He I, you have to deal with the presence of the
   nearby Si I line because the He I multiplet is in the wings of this photospheric line.
   In this version of the code you can use SIR to simuultaneously invert the photospheric
   and chromospheric lines, so you do not have to worry for decontaminating from the
   Si I line. You can even do that with the telluric contamination nearby.		
   If the data is off-limb, things are typically easier because there is no continuum but		
   sometimes there is some stray-light that can give you a headache. In this case, 		
   the input should be normalized by the peak emission.		
		
-  **Wavelength calibration**: the data has to be wavelength calibrated. How to do it depends 		
   on whether you want an absolute calibration of velocities or not. If you want such absolute scale, the best is to do the wavelength		
   calibration using telluric lines and then transform everything to the Sun using the relative		
   velocity between the observed region and the Earth. If not, maybe using some weak 		
   surrounding photospheric lines is enough. 
		
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
   geometry in mind and try to adapt it to your observations. See :ref:`equations` for more information.
		
