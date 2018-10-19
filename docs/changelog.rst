Changelog
=========

Version 2018.9.19
-----------------
 - No external files needed for SIR
 - Coupling of parameters
 - Use other optimization algorithms
 - Check that initial parameters are inside the range if not inverted
 - Masks
 - Save chi2 for each cycle
 - Comment on virtual environments
 - Explain better how atmospheres are defined in Hazel
 - Randomization
 - Error bars
 - Allow to invert magnetic fields in the LOS or vertical reference frames
 - Instrumental profile
 - More tests added
 - Output AIC and BIC


To do
-----

 - Add even more tests
 - Sparse regularization of maps : change order of loops -> 1 LM (x,y) -> regularization
 - Regularization of parameters
 - In synthesis, all atmospheres have to have the same number of pixels
 - Add utilities for generating new SIR models
 - Show an example for inverting a random distribution of points to be inverted and then start from an interpolation of these maps to be closer to the solution
 - Develop tools to easily manage input/output models so that multiinversions are easy to do.
 - Check that output and input models are compatible so that reinversions can be easily done
 - Precalculate coefficients for the emission/absorption coefficients using sparse matrices and use them later on the calculation. This will avoid many of the loops that are, by construction, zero.
 - Programmatic inversions