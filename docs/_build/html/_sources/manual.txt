Introduction
============

Description
-----------

 (an acronym for HAnle and ZEeman Light) is a computer program for the
synthesis and inversion of Stokes profiles caused by the joint action of
atomic level polarization and the Hanle and Zeeman effects. It is based
on the quantum theory of spectral line polarization, which takes into
account rigorously all the relevant physical mechanisms and ingredients:
optical pumping, atomic level polarization, level crossings and
repulsions, Zeeman, Paschen-Back and Hanle effects. The code is written
in standard Fortran 90. Its parameters are passed using four
configuration files that can be manually edited. These configuration
files are heavily commented, so that their edition should be an easy
task. In any case, two front-ends coded in IDL are given as a part of
the distribution in order to facilitate a user-friendly execution of the
program. A parallel version of the code using Message Passing Interface
(MPI) is also available. This manual considers both distributions.

Credits
-------

The code has grown since the first version thanks to the suggestions of
many people. We thank Rebecca Centeno Elliot, Yukio Katsukawa, Marian
Martínez González, Rafael Manso Sainz and Tom Schad for their help on
testing the code and proposing (and partially coding, in some cases)
some of the options of the code.

Uncompressing and compiling 
============================

Serial version
--------------

The package comes in a single compressed file ``hazel.tar.gz``. After
unpacking with ``tar zxvf hazel.tar.gz``, the  directory will contain
the following subdirectories:

#. Source contains the Fortran 90 sources and a makefile that can be
   used to build the binary file.

#. Run contains a directory tree structure with the appropriate
   configuration files to run the code in command line mode.

#. Widget\_Synth contains all the files that are needed to run the IDL
   front-end for the synthesis problem.

#. Widget\_Inv contains all the files that are needed to run the IDL
   front-end for the inversion problem.

#. IDL\_routines contains some IDL routines that are needed by the
   front-ends.

#. Manual contains this manual.

The code has been tested on Linux platforms using the Intel Fortran
Compiler (``ifort``) and the free GFortran compiler. The source code is
in the ``Source/`` directory. The compilation is performed with the
supplied ``makefile``. It is quite simple and easy to modify, and
contains additional comments about compiling. The default compiler is
the ``ifort``, although you can use any other compiler through the
variable ``COMPILER``. In order to obtain the executable file, just
type:

::

           make all

After compiling and linking, the executable is copied to the  ``Run/``,
``Widget_Synth/`` and ``Widget_Inv/`` directories. Running the program
in the ``Run/`` directory should produce the correct output depending on
the exact form of the input files.

The generated object and module files can be cleaned typing:

::

           make clean

Parallel version
----------------

The package also decompresses the  directory tree that will contain the
following subdirectories:

#. SourceMPI contains the Fortran 90 sources and a makefile that can be
   used to build the binary file.

#. RunMPI contains a directory tree structure with the appropriate
   configuration files to run the code in command line mode.

The source code is in the ``SourceMPI/`` directory. The compilation
depends on the precompiled library NetCDF [2]_ for reading and writing
output files. NetCDF is a standard for platform independent binary files
that you need to have installed in your system. The compilation is
performed with the supplied ``makefile``. It is quite simple and easy to
modify, and contains additional comments about compiling. The default
compiler is ``mpif90``, although you can use any other compiler through
the variable ``COMPILER``. The variables ``NETCDF_INCLUDE`` and
``NETCDF_LIB`` have to point to the ``include`` and ``lib`` directories
of the NetCDF distribution.

The code makes use of the MPI package for parallelization, so it has to
be installed on your system. In order to obtain the executable file (for
instance for the Intel compiler), just type:

::

           make -f makefile.Intel

Modify the ``makefile`` to point the variables to the correct libraries
and include files. After compiling and linking, the executable is copied
to the  ``RunMPI/`` directory, where the code is run. Running the
program in the ``RunMPI/`` directory should produce the correct output
depending on the exact form of the input files.

The generated object and module files can be cleaned typing:

::

           make clean

The code is run from the ``RunMPI`` directory. Use your MPI launcher to
select the number of processors. For example:

::

    mpiexec -n 50 hazel_mpi config_inversion.dat 2000 5000

The code admits up to three command line parameters:

-  Filename with the main configuration file.

-  Starting pixel of the inversion. This is used if you want to rerun
   the inversion of some pixels.

-  Final pixel of the inversion. This is used if you want to rerun the
   inversion of some pixels.

See §[sec:phazel\_files] for details on the input files.

New input file
==============

In previous versions of the code, the code was controlled with four
configuration files. The updated version of the code is now controlled
by one human-readable configuration file. The code can still be run
using the old configuration files but this option will be discontinued
in the future so it is advisable to to the shift to this configuration
file. In order to use this option, you need to have the ``configparser``
package installed in your system. It can be downloaded from
``http://www.voidspace.org.uk/python/configobj.html``. The serial code
is now run using:

::

    ./run.py conf.ini

and the parallel code is run with

::

    ./run.py conf.ini nProcessors

An example of the file, that is self-explanatory, is:

::

    # Hazel configuration File

    #####################
    # General information
    #####################

    [Files]
    Input model file = 'ATOMS/helium.mod'
    File with observations = 'OBSERVATION/test_2comp.prof'
    File with inverted profiles = 'test.inversion'
    File with inverted parameters = 'test.parameters'

    [Working mode]
    Action = 'inversion'                        # 'synthesis' or 'inversion'
    Verbose = no
    Linear system solver = 'LU'                 # 'LU' or 'CG'
    Stopping volume for DIRECT = 0.001

    [General parameters]
    Synthesis mode = 'exact'                    # 'thin' or 'exact'
    Include stimulated emission = yes
    Include magnetic field  = yes
    Include Paschen-Back effect = yes
    Include atomic level polarization = yes
    Include magneto-optical effects in the RT = yes
    Include stimulated emission in the RT = yes
    Multiplet = 10830                           # 10830, 5876, 7065, 3889 A
    Line-of-sight angles = 0.0, 0.0, 90.0       # theta, chi, gamma deg
    Wavelength axis = -3.0, 2.5, 200            # Minimum, maximum and number of grid points

    #####################
    # Synthesis parameters
    #####################
    [Synthesis]
    Number of slabs = '1'                       # '1' -> single slab, '1+1' -> two slabs with same field, '1+1B' -> 2 slabs with different field, '2' -> two slabs added with a filling factor
    Boundary condition = 4.098e-5, 0.0, 0.0, 0.0      # I0, Q0, U0, V0
    a = 0.0
    height = 3.0                                # Real height if positive, apparent height if negative arcsec
    ff = 0.0
        [[Slab 1]]
        B =         0.0         # G
        thetaB =    0.0         # deg
        chiB =      0.0         # deg
        vdopp =     8.0         # km/s
        tau =       1.0
        vmac =      0.0         # Positive is redshift km/s
        beta =      1.0
        [[Slab 2]]
        B =         0.0         # G
        thetaB =    0.0         # deg
        chiB =      0.0         # deg
        vdopp =     0.0         # km/s
        tau =       0.0
        vmac =      0.0         # Positive is redshift km/s
        beta =      1.0

    #####################
    # Ranges for the DIRECT method [min, max]
    #####################
    [Ranges]
    a =             0,0.5
    ff =            0.0,1.0
        [[Slab 1]]
        B =         800,1100
        thetaB =    0,180
        chiB =      0,180
        vdopp =     2,12
        tau =       0.1,2
        vmac =      -5,5
        beta =      0.5,2
        [[Slab 2]]
        B =         800,1100
        thetaB =    0,180
        chiB =      0,180
        vdopp =     2,12
        tau =       0.1,2
        vmac =      -5,5
        beta =      0.5,2
        
    #####################
    # Parameters to invert
    #####################
    [Inversion]
    Iterations in LM = 20
    Number of cycles = 4
    Inversion modes = 'DIRECT', 'LM', 'DIRECT', 'LM'        # 'DIRECT' for DIRECT algorithm and 'LM' for Levenberg-Marquardt
        [[Cycles]]
        a =             1, 1, 0, 0
        ff =            0, 0, 0, 0
            [[[Slab 1]]]
            B =         0, 0, 1, 1
            thetaB =    0, 0, 1, 1
            chiB =      0, 0, 1, 1
            vdopp =     1, 1, 0, 0
            tau =       1, 1, 0, 0
            vmac =      1, 1, 0, 0
            beta =      0, 0, 0, 0
            [[[Slab 2]]]
            B =         0, 0, 0, 0
            thetaB =    0, 0, 0, 0
            chiB =      0, 0, 0, 0
            vdopp =     0, 0, 0, 0
            tau =       0, 0, 0, 0
            vmac =      0, 0, 0, 0
            beta =      0, 0, 0, 0
        [[Weights]]
            Stokes I =  1.0, 1.0, 1.0, 1.0
            Stokes Q =  0.0, 0.0, 1.0, 1.0
            Stokes U =  0.0, 0.0, 1.0, 1.0
            Stokes V =  0.0, 0.0, 1.0, 1.0

Input files
===========

 is controlled via four configuration files. All configuration files are
fully commented, so that changing any parameter should be an easy task.
In the following, we describe them step by step.

``config_inversion.dat``
------------------------

This file can be considered as the main configuration file and it is the
only one that has to have a fixed name. This file is used to indicate
the names of the input files, the names of the output files, verbosity
level and to decide whether  is to be applied to work in synthesis or
inversion mode. Using the example included in the present version of ,
we analyze one by one all the inputs.

::

    # Input model file
    'ATOMS/helium.mod'

Definition of the file with the atomic model. See §[sec:atomic\_model]
for an explanation of the file format.

::

    # Initial parameters file
    'init_parameters.dat'

Definition of the file with the initial parameters of the problem. The
values of the parameters in this file are taken as initial values for
the inversion or for the synthesis. See §[sec:init\_parameters] for a
detailed description of the file.

::

    # Range of parameters for the DIRECT method
    'direct_range.dat'

This file is used to define the lower and upper limits of the intervals
inside which the DIRECT method searches for the minimum of the
:math:`\chi^2` function. See §[sec:direct\_range] for details.

::

    # Output for the upper level rho^K_Q(J,J') in the vertical reference frame
    'ATOMIC_POL/vertical_upper.rho'

    # Output for the lower level rho^K_Q(J,J') in the vertical reference frame
    'ATOMIC_POL/vertical_lower.rho'

    # Output for the upper level rho^K_Q(J,J') in the mag. field reference frame
    'ATOMIC_POL/magnetic_upper.rho'

    # Output for the lower level rho^K_Q(J,J') in the mag. field reference frame
    'ATOMIC_POL/magnetic_lower.rho'

The previous lines define the output files where the spherical tensor
components of the density matrix are saved. Note that the code stores
only the density matrix elements of the upper and lower level of the
desired transition. The elements of the atomic density matrix depend on
the chosen reference system, and the two most desired reference systems
are the one in which the quantization axis is chosen along the solar
local vertical direction and the one in which the quantization axis is
chosen along the magnetic field vector.

::

    # Output absorption/emission coefficients
    'INVERTED/rtcoef.emer'

    # Output absorption/emission coefficients neglecting atomic polarization
    'INVERTED/rtcoef_noatompol.emer'

The emission coefficients :math:`\epsilon_{I,Q,U,V}`, the absorption
coefficients :math:`\eta_{I,Q,U,V}` and the anomalous dispersion
coefficients :math:`\rho_{Q,U,V}` for each wavelength point are saved in
these files. The first file includes the effects of atomic level
polarization, while the second one neglects its influence.

::

    # File with the observed profiles
    'OBSERVATION/test.prof'

When using the code in the inversion mode, this file is the one used for
the input of the observed Stokes profiles. The format of this file
depends on which version of the code is used. For , it is very simple.
The first line indicates the number of wavelength points and the
normalization (use ’cont’ or ’peak’). Then, a table with nine columns
gives the value of the wavelength shift with respect to the center of
the multiplet, the Stokes vector at each wavelength normalized to the
maximum intensity, and an estimation of the noise standard deviation at
each wavelength normalized to the maximum intensity. See the example
file contained in the  distribution for more details. Note that these
lines have to be present in the input file even if  is used in synthesis
mode.

When using , the input file is more complicated and is described in
§[sec:phazel\_files].

::

    # File with the inverted profiles
    'test.inversion'

    # File with the parameters from the inversion
    'test.parameters'

The final Stokes profiles resulting from the synthesis or inversion
options is saved in the file indicated in the first line. The format is
the same as that explained for the file containing the observation. When
 is run in inversion mode, the final inferred parameters of the model
are saved in the file indicated in the second line. Again, for  the
output files are described in §[sec:phazel\_files].

::

    # File that sets the parameters to invert
    'invert_parameters.dat'

This file defines which parameters to invert in the inversion mode,
together with the algorithm to be used in each cycle and the weight used
for each Stokes parameter.

::

    # Verbose mode (0-> no, 1-> yes)
    0

Flag to connect or disconnect the verbose mode. For the inversion of
Stokes profiles affected by atomic level polarization it is sometimes
useful to turn the verbose mode on for analyzing the process of the code
while calculating.

::

    # Linear system solver (0-> LU, 1-> CG)
    0

This flag is used to choose the algorithm that solves the linear system
of statistical equilibrium equations. For relatively simple models, the
LU decomposition does a very good job in terms of speed. If the number
of unknowns (i.e., of :math:`\rho^K_Q(J,J')` elements) turns out to be
of the order of or larger than :math:`10^3`, conjugate gradients (CG)
methods are a much better option. We recommend to use the LU
decomposition when possible and move to the CG solution only when
necessary. The CG solution are based on routines developed by Dr. Mark
K. Seager from Lawrence Livermore National Lab.

::

    # Optically thin (0), slab no-MO (1), M-E (2), slab DELOPAR (3), 
                           simplified slab (4), exact slab (5)
    5

This flag is used to choose the level of approximation for the solution
of the radiative transfer equation. The meaning of each option is
explained below in §[sec:radiative\_transfer].

::

    # Synthesis mode -> 0 , Inversion mode -> 1
    0

This flag controls the working mode of the code (synthesis or
inversion).

``init_parameters.dat``
-----------------------

This important file establishes the parameters of the model, together
with the definition of the scattering geometry. It includes also flags
to turn on or discard different physical mechanisms. In the synthesis
mode, the values in this file are used to carry out the synthesis. In
the inversion mode, the values in this file are chosen as initial
conditions for the inversion for those parameters that are left free.
For those that are left fixed, the code uses the values defined in this
file. We explain them step by step.

::

    # Include stimulated emission (0-> no, 1-> yes)
    1

This flag is used to take into account or discard the effect of
stimulated emission in the emergent Stokes profiles. Although stimulated
emission is negligible for most solar it can be of importance for very
strong radiation fields. We recommend to use always 1 since the
computational time is barely affected by this flag.

::

    # Include magnetic field (0-> no, 1-> yes)
    1

This flag is used to slightly reduce the computational work for the
non-magnetic case because, if set to zero, the magnetic kernel [see Eq.
([eq:see])] is not calculated.

::

    # Include depolarization rates (0-> no, 1-> yes)
    0

    # Value of delta if depol. rates are included (not used if prev. value = 0)
    1.d14

In the present version of  it is possible to include the effect of
depolarizing collisions only in the ground level of the atomic system.
In case the effect of collisions is to be accounted for, set the first
parameter to 1 and give the collisional rate in the next parameter in
units of s\ :math:`^{-1}`.

::

    # Include Paschen-Back effect (0-> no, 1-> yes)
    1

The effect of a magnetic field on the energy levels of the atomic system
can be calculated under the approximation of the linear Zeeman effect or
in the general case of the intermediate Paschen-Back effect. If this
flag is set to 0, the approximation of the linear Zeeman effect is used
and no perturbations between different :math:`J` levels of a term are
taken into account. If the flag is set to 1, the general theory of the
Paschen-Back effect is used to calculate the wavelength positions and
the strengths of the :math:`\pi` and :math:`\sigma` components. The
difference in the computational work between both approaches is rather
small.

::

    # Number of slabs (1-> 1 slab, 2-> 2 slabs with same B, 
    3-> 2 slabs with different B, -2 -> 2 slabs with filling factor)

 can be used using one slab (option 1) of constant physical properties
or two (options 2 and 3 and -2). The difference between options 2 and 3
is that option 2 considers both slabs to have exactly the same field
while option 3 considers two different fields. As a consequence, the
computing time is smaller in option 2. In both options, the second slab
is placed in front of the first one, so that the boundary condition of
the second slab is the emergent radiation from the first. In option -2,
the radiation emerging from both slabs is added weighted with a filling
factor, which is indicated below.

::

    # Magnetic field strength [G], thetaB [degrees], chiB [degrees]
    0.3d0 90.d0 90.d0

The magnetic field vector is defined here. The strength in G and the
inclination and azimuth angles in degrees define the magnetic field
vector. The angles are defined with respect to the vertical direction in
the atmosphere, as shown in Fig. [fig:geometry]. Note that, if the
azimuth of the field is set to 999, the random azimuth solution is
obtained following the strategy explained in Appendix C of
:raw-latex:`\cite{belluzzi07}`. If two slabs are used (setting option 3
or -2 above), put the two field vectors next to each one in the format
:math:`(B,\theta_B,\chi_B)_1 (B,\theta_B,\chi_B)_2`.

::

    # Apparent height (if <0) or real height (if >0) of the atoms in arcsec
    3.d0

The tensors :math:`J^0_0` and :math:`J^2_0` that quantify the mean
intensity of the radiation field and its anisotropy are calculated
assuming a standard solar center-to-limb variation (CLV) and taking into
account geometrical effects. This parameter gives the height at which
the slab of atoms is placed with respect to the surface of the Sun.

::

    # Optical depth of the slab in the maximum of I (slab) or strength of the line (ME)
    1.0d0

This quantity is the optical depth of the slab at the wavelength
position of the maximum absorption or emission in Stokes :math:`I`. For
example, for the 10830 Å multiplet of He i, this is the position of the
red blended component. If two slabs with option 2 or 3 are used, put the
two optical depths together. If option -2 is used, then add the filling
factor as a third number.

::

    # Source function increase
    1.d0

The source function of the slab will be multiplied by this number. This
is a way to generate lines in emission even when the slab is seen on the
solar disk. If two components (one after the other) are used, this
number only modifies the source function of the second component. This
allows us to simulate self-absorption in the code.

::

    # Boundary Stokes parameters (I0,Q0,U0,V0)
    4.098d-5 0.d0 0.d0 0.d0

Boundary conditions for the Stokes vector used in the solution of the
radiative transfer equation. If the radiation field is the photospheric
continuum, the IDL routine ``IDL_routines/solar_field.pro`` can be used
to return an estimation.

::

    # Transition where to compute the emergent Stokes profiles
    1

From the transitions defined in the atomic model, the code calculates
the emergent Stokes profiles for the chosen transition. For the moment,
only one transition at a time is allowed. We plan to extend this to
synthesize several lines.

::

    # Include atomic level polarization? (0-> no, 1-> yes)
    1

The synthesis or inversion options can be used taking into account or
neglecting the presence of atomic level polarization. This flag controls
it.

::

    # Observation angle with respect to the local solar vertical theta,chi,gamma [degrees]
    0.d0 0.d0 90.d0

The line-of-sight direction is defined using the angles described in
Fig. [fig:geometry]. All angles are given in degrees.

::

    # Wavelength axis: minimum, maximum and number of grid points
    -3.d0 2.5d0 200

In case the code is run in synthesis mode, this line is used to set the
lower and upper limits (in cm\ :math:`^{-1}`) of the wavelength axis.
The last parameter gives the number of wavelength points to be used. In
the inversion mode, the wavelength axis is chosen automatically from the
observation and these numbers are overridden.

::

    # Line wavelength [A], Doppler velocity [km/s] and damping [a]
    10829.0911d0   6.5d0   0.d0

This line is used to define the wavelength of the multiplet (wavelength
of the :math:`(L,S) \to (L',S')` transition), the Doppler width of the
line in km s\ :math:`^{-1}` and the reduced damping constant. If two
slabs (through options 3 or -2) are used, add the Doppler width of the
second component next to the first one. Concerning the reduced damping
constant, if its value is negative, it is computed using the natural
damping and using the Doppler broadening. The absolute value of the
input value is used then as an enhancement factor (so you should use
:math:`-1` is you want to use the natural width).

::

    # Macroscopic velocity [km/s] (>0 is a redshift)
    0.d0

This defines the wavelength shift produced by the presence of a bulk
motion of the plasma. Note that positive velocities imply redshifts. If
two components (options 2, 3 or -2) are used, put the two bulk
velocities.

::

    # Include magneto-optical effects in the RT
    1

It is possible to include (1) or neglect (0) the influence of the
anomalous dispersion coefficients :math:`\rho_{Q,U,V}` in the
calculation of the emergent Stokes profiles.

::

    # Include stimulated emission in the RT
    1

This flag controls whether we include (1) or neglect (0) the influence
of the stimulated emission in the calculation of the emergent Stokes
profiles.

``direct_range.dat``
--------------------

The DIRECT global optimization method is used to give a first estimation
of the parameters from which the Levenberg-Marquardt method is applied
to locate the minimum of the :math:`\chi^2` surface. The behavior of the
DIRECT method is controlled with this file, in which we must specify the
upper and lower limits of the model parameters, together with details
about the stopping criterion. In the following, we describe all the
options in detail.

::

    # Output file
    'direct.location'

The DIRECT method tries to evaluate the merit function :math:`\chi^2` as
few times as possible. The code saves in this file the values of the
parameters at which the algorithm has carried out the evaluation of the
merit function. This can be useful for analyzing the presence of
ambiguities. In this case, the method will clearly mark the position of
the possible solutions by evaluating the merit function more times in
the surroundings of the compatible solutions. Note that this lines are
absent on the  configuration file.

::

    # Maximum number of function evaluations (<0 -> don't use this criteria)
    -1

    # Reduction in the volume (<0 -> don't use this criteria, typically 0.01)
    0.001

The previous two lines are used to indicate the stopping criterion for
the DIRECT method. An early stop will probably give a first estimation
of the solution that is far from the final result. Letting the code run
for many iterations may degrade too much the computing time because of
the poor local convergence properties of the DIRECT scheme. The first
option permits the user to stop after a fixed number of evaluations of
the merit function. The second option permits the user to stop when the
ratio between the hypervolume where the global minimum is located and
the original hypervolume is smaller than the given threshold. We have
verified that 0.001 gives very good results. Setting one of the two
parameters to values :math:`< 0` will disconnect it.

::

    # Magnetic field (0-Bmax)
    800.d0  1100.d0

    # thetab  (0 .. 180)
    30.d0  180.d0

    # chib (0 .. 180)
    -180.d0  0.d0

    # vdopp (0 .. 20)
    2.d0  7.d0

    # dtau (0 .. 5)
    0.d0  1.d0

    # delta_collision (0 .. 18)
    0.d0  18.d0

    # vmacro (-10 .. 10)
    -10.d0  10.d0

    # damping (0 .. 4)
    0.d0  4.d0

    # beta (0 .. 10)
    0.d0  1.d0

    # height (0 .. 100)
    0.d0  100.d0

    # dtau2 (0 .. 5)
    0.d0  2.d0

    # vmacro2 (-10 .. 10)
    25.d0  35.d0

    # Magnetic field 2 (0-Bmax)
    800.d0  1100.d0

    # thetab 2 (0 .. 180)
    30.d0  180.d0

    # chib 2 (0 .. 180)
    -180.d0  0.d0

    # vdopp 2 (0 .. 20)
    2.d0  12.d0

The previous lines define the space of parameters where the DIRECT
method will look for the global minimum.

``invert_parameters.dat``
-------------------------

This file is used to set the behavior of the inversion mode: the
structure of the inversion cycle, setting the free and the fixed
parameters.

::

    # Maximum number of iterations
    20

This parameter sets the maximum number of Levenberg-Marquardt (LM)
iterations to be carried out in each cycle. Sometimes the LM scheme
stops before reaching the maximum number of iterations because the
relative change in the parameters from one iteration to the next is
below 10\ :math:`^{-4}`.

::

    # Number of cycles
    2

The optimal iteration scheme is composed of combinations of cycles. In
the first cycle, the DIRECT method is used to give a first estimation of
the solution. In the second cycle, the LM method is used to refine the
solution until arriving to the final one. This parameter sets the number
of cycles used.

::

    # Invert the magnetic field strength
    1 1 1 1

    # Invert the magnetic field inclination
    1 1 1 1

    # Invert the magnetic field azimuth
    1 1 0 0

    # Invert the Doppler width
    0 0 0 0

    # Invert the optical depth or strength of the line
    0 0 0 0

    # Invert the D^2 of the lower level
    0 0 0 0

    # Invert the macroscopic velocity
    0 0 0 0

    # Invert the damping
    0 0 0 0

    # Invert the source function gradient
    0 0 0 0

    # Invert the height of the He atoms
    0 0 0 0

    # Invert the optical depth or strength of the line of component 2
    0 0 0 0

    # Invert the macroscopic velocity of component 2
    0 0 0 0

    # Invert the magnetic field strength of component 2
    0 0 1 1

    # Invert the magnetic field inclination of component 2
    0 0 1 1

    # Invert the magnetic field azimuth of component 2
    0 0 1 1

    # Invert the Doppler width of component 2
    0 0 0 0

Depending on the number of cycles, the previous lines define whether a
parameter is inverted (setting a 1 in the corresponding cycle) or kept
fixed to the value given in the ``init_parameters.dat`` file (setting a
0 in the corresponding cycle). The number of 0s/1s in each line has to
be larger or equal to the number of cycles.

::

    # Weights for Stokes I in each cycle
    1.d0 1.d0 1.d0 1.d0

    # Weights for Stokes Q in each cycle
    1.d0 1.d0 1.d0 1.d0

    # Weights for Stokes U in each cycle
    1.d0 1.d0 1.d0 1.d0

    # Weights for Stokes V in each cycle
    1.d0 1.d0 1.d0 1.d0

Since the inversion is based on the gradient descent minimization of the
:math:`\chi^2` merit function and not on sampling methods, it is
important to modify sometimes the weight of each Stokes vector in order
to increase the sensitivity of the :math:`\chi^2`-function to some model
parameters. The code allows to change the relative weight of each Stokes
vector in each cycle.

::

    # Inversion modes (1-> LM, 2-> DIRECT, 3-> PIKAIA)
    2 1 2 1

The optimization method used in each cycle is set in this line. Note
that the scheme DIRECT+LM has been empirically proved to be quite
optimal. The possibility to use genetic optimization based on the Pikaia
algorithm is still in a preliminary phase. However, the large number of
function evaluations that any genetic algorithm needs makes it difficult
to beat the DIRECT+LM combination.






Ambiguities
-----------

You have to remember that the results of Hazel are potentially affected
by ambiguities and you have to take them into account. There is an
utility written in IDL that, given an inverted map, obtains all the
other solutions which are ambiguous in the saturation regime. This can
be called, including the appropriate paths and discarding the final
``.nc`` extension, by:

::

    IDL> disamb, 'file_with_inversions', 'file_with_observations', angleObs

where ``angleObs`` is the observation angle :math:`\theta` (so that it
is 90\ :math:`^\circ` for an observation exactly at the limb. This
program can be called with the additional ``/gen_files_inversion``,
which then generates a set of observations, configuration files and a
file to run . This is useful in case the line is not in the saturation
regime. In this case, the ambiguous solutions that are found by the code
are not strictly valid and one should refine them with a final LM cycle
in which :math:`B`, :math:`\theta_B` and :math:`\chi_B` are left free.
The solution to the ambiguities in the saturation regime is shown in
Section [sec:ambiguities].

Calling Hazel from Python
=========================

We have developed a wrapper to allow the user to call the synthesis
routines of Hazel in Python. To do so, just enter into the directory
``SourcePy`` and type

::

    python setup.py build_ext --inplace

and a library ``pyhazel.so`` will be generated (and also copied to the
directory ``RunPy``. In this very same directory you can see the
``test.py`` file that shows how to call the code to wrapper.

Basic Equations
===============

We consider a constant-property slab of atoms, located at a height
:math:`h` above the visible solar “surface", in the presence of a
deterministic magnetic field of arbitrary strength :math:`B`,
inclination :math:`\theta_B` and azimuth :math:`\chi_B` (see Fig. 1).
The slab’s optical thickness at the wavelength and line of sight under
consideration is :math:`\tau`. We assume that all the atoms inside this
slab are illuminated from below by the photospheric solar continuum
radiation field, whose center-to-limb variation has been tabulated by
:raw-latex:`\cite{pierce00}`. The ensuing anisotropic radiation pumping
produces population imbalances and quantum coherences between pairs of
magnetic sublevels, even among those pertaining to the different
:math:`J`-levels of the adopted atomic model. This atomic level
polarization and the Zeeman-induced wavelength shifts between the
:math:`\pi` (:math:`\Delta{M}=M_u-M_l=0`), :math:`\sigma_{\rm blue}`
(:math:`\Delta{M}=+1`) and :math:`\sigma_{\rm red}`
(:math:`\Delta{M}=-1`) transitions produce polarization in the emergent
spectral line radiation.

In order to facilitate the understanding of the code, in the following
we summarize the basic equations which allow us to calculate the
spectral line polarization taking rigorously into account the joint
action of atomic level polarization and the Hanle and Zeeman effects. To
this end, we have applied the quantum theory of spectral line
polarization, which is described in great detail in the monograph by
:raw-latex:`\cite{landi_landolfi04}`. We have also applied several
methods of solution of the Stokes-vector transfer equation, some of
which can be considered as particular cases of the two general methods
explained in §6 of :raw-latex:`\cite{trujillo03}`.

.. figure:: f1.eps
   :alt: The geometry for the scattering event. The :math:`Z`-axis is
   placed along the vertical to the solar atmosphere. The magnetic field
   vector, :math:`\mathbf{B}`, is characterized by its modulus
   :math:`B`, the inclination angle :math:`\theta_B` and the azimuth
   :math:`\chi_B`. The line-of-sight, indicated by the unit vector
   :math:`\mathbf{\Omega}`, is characterized by the two angles
   :math:`\theta` and :math:`\chi`. The reference direction for Stokes
   :math:`Q` is defined by the vector :math:`\mathbf{e}_1` on the plane
   perpendicular to the line-of-sight. This vector makes an angle
   :math:`\gamma` with respect to the plane formed by the vertical and
   the line-of-sight. In the figures showing examples of the emergent
   Stokes profiles, our choice for the positive reference direction for
   Stokes :math:`Q` is :math:`\gamma=90^\circ`, unless otherwise stated.
   For off-limb observations, we have :math:`\theta=90^\circ`, while for
   observations on the solar disk, we have :math:`\theta<90^\circ`. Note
   also that :math:`\chi` is generally taken to be :math:`0^\circ`.
   [fig:geometry]

   The geometry for the scattering event. The :math:`Z`-axis is placed
   along the vertical to the solar atmosphere. The magnetic field
   vector, :math:`\mathbf{B}`, is characterized by its modulus
   :math:`B`, the inclination angle :math:`\theta_B` and the azimuth
   :math:`\chi_B`. The line-of-sight, indicated by the unit vector
   :math:`\mathbf{\Omega}`, is characterized by the two angles
   :math:`\theta` and :math:`\chi`. The reference direction for Stokes
   :math:`Q` is defined by the vector :math:`\mathbf{e}_1` on the plane
   perpendicular to the line-of-sight. This vector makes an angle
   :math:`\gamma` with respect to the plane formed by the vertical and
   the line-of-sight. In the figures showing examples of the emergent
   Stokes profiles, our choice for the positive reference direction for
   Stokes :math:`Q` is :math:`\gamma=90^\circ`, unless otherwise stated.
   For off-limb observations, we have :math:`\theta=90^\circ`, while for
   observations on the solar disk, we have :math:`\theta<90^\circ`. Note
   also that :math:`\chi` is generally taken to be :math:`0^\circ`.
   [fig:geometry]

The radiative transfer approach
-------------------------------

The emergent Stokes vector
:math:`\mathbf{I}(\nu,\mathbf{\Omega})=(I,Q,U,V)^{\dag}` (with
:math:`\dag`\ =transpose, :math:`\nu` the frequency and
:math:`\mathbf{\Omega}` the unit vector indicating the direction of
propagation of the ray) is obtained by solving the radiative transfer
equation

.. math::

   \frac{d}{ds}\mathbf{I}(\nu,\mathbf{\Omega}) =
   \bm{\epsilon}(\nu,\mathbf{\Omega}) - \mathbf{K}(\nu,\mathbf{\Omega}) 
   \mathbf{I}(\nu,\mathbf{\Omega}),
   \label{eq:rad_transfer}

where :math:`s` is the geometrical distance along the ray under
consideration,
:math:`\bm{\epsilon}(\nu,\mathbf{\Omega})=({\epsilon}_I,{\epsilon}_Q,{\epsilon
}_U,{\epsilon}_V)^{\dag}` is the emission vector and

.. math::

   \mathbf{K} = \left( \begin{array}{cccc}
   \eta_I & \eta_Q & \eta_U & \eta_V \\
   \eta_Q & \eta_I & \rho_V & -\rho_U \\
   \eta_U & -\rho_V & \eta_I & \rho_Q \\
   \eta_V & \rho_U & -\rho_Q & \eta_I
   \end{array} \right)
   \label{eq:propagation}

 is the propagation matrix. Alternatively, introducing the optical
distance along the ray, :math:`{\rm d}{\tau}=-{\eta_I}{\rm d}s`, one can
write the Stokes-vector transfer Eq. ([eq:rad\_transfer]) in the
following two ways:

-  The first one, whose formal solution requires the use of the
   evolution operator introduced by :raw-latex:`\cite{landi_landi85}`,
   is

   .. math::

      {{d}\over{d{\tau}}}{\bf I}\,=\,{\bf K}^{*}
      {\bf I}\,-\,{\bf S}, 
      \label{eq:rad_transfer_peo}

    where :math:`{\bf K}^{*}={\bf K}/{\eta_I}` and
   :math:`{\bf S}=\bm{\epsilon}/{\eta_I}`. The formal solution of this
   equation can be seen in eq. (23) of :raw-latex:`\cite{trujillo03}`.

-  The second one, whose formal solution does not require the use of the
   above-mentioned evolution operator is
   :raw-latex:`\citep[e.g.,][]{rees_delo89}`

   .. math::

      {{d}\over{d{\tau}}}{\bf I}\,=\,{\bf I}\,-\,{\bf S}_{\rm eff},  
      \label{eq:rad_transfer_delo}

    where the effective source-function vector
   :math:`\,{\bf S}_{\rm eff}\,=\,{\bf S}\,-\,
   {\bf K}^{'}{\bf I},\,\,\,` being
   :math:`\,{\bf K}^{'}={\bf K}^{*}-{\bf 1}` (with :math:`\bf 1` the
   unit matrix). The formal solution of this equation can be seen in eq.
   (26) of :raw-latex:`\cite{trujillo03}`.

Once the coefficients :math:`\epsilon_I` and :math:`\epsilon_X` (with
:math:`X=Q,U,V`) of the emission vector and the coefficients
:math:`\eta_I`, :math:`\eta_X`, and :math:`\rho_X` of the
:math:`4\times4` propagation matrix are known at each point within the
medium it is possible to solve formally Eq. ([eq:rad\_transfer\_peo]) or
Eq. ([eq:rad\_transfer\_delo]) for obtaining the emergent Stokes
profiles for any desired line of sight. Our computer program considers
the following levels of sophistication for the solution of the radiative
transfer equation:

-  *Numerical Solutions*. The most general case, where the properties of
   the slab vary along the ray path, has to be solved numerically. To
   this end, two efficient and accurate methods of solution of the
   Stokes-vector transfer equation are those proposed by
   :raw-latex:`\cite{trujillo03}` (see his eqs. (24) and (27),
   respectively). The starting points for the development of these two
   numerical methods were Eq. ([eq:rad\_transfer\_peo]) and Eq.
   ([eq:rad\_transfer\_delo]), respectively. Both methods can be
   considered as generalizations, to the Stokes-vector transfer case, of
   the well-known short characteristics method for the solution of the
   standard (scalar) transfer equation.

-  *Exact analytical solution of the problem of a constant-property slab
   including the magneto-optical terms of the propagation matrix*. For
   the general case of a constant-property slab of arbitrary optical
   thickness we actually have the following analytical solution, which
   can be easily obtained as a particular case of eq. (24) of
   :raw-latex:`\cite{trujillo03}`:

   .. math::

      {\bf I}={\rm e}^{-{\mathbf{K}^{*}}\tau}\,{\bf I}_{\rm sun}\,+\,\left[{\mathbf{K}^{*}}\right]^{-1}\,
      \left( \mathbf{1} - {\rm e}^{-{\mathbf{K}^{*}}\tau} \right) \,\mathbf{S},
      \label{eq:slab_peo}

   where :math:`\mathbf{I}_{\rm sun}` is the Stokes vector that
   illuminates the slab’s boundary that is most distant from the
   observer. We point out that the exponential of the propagation matrix
   :math:`{\mathbf{K}^{*}}` has an analytical expression similar to eq.
   (8.23) in :raw-latex:`\cite{landi_landolfi04}`.

-  *Approximate analytical solution of the problem of a
   constant-property slab including the magneto-optical terms of the
   propagation matrix*. An approximate analytical solution to the
   constant-property slab problem can be easily obtained as a particular
   case of eq. (27) of :raw-latex:`\cite{trujillo03}`:

   .. math::

      \mathbf{I} = \left[ \mathbf{1}+\Psi_0 \mathbf{K}' \right]^{-1} \left[ \left(
      e^{-\tau} \mathbf{1} - \Psi_M \mathbf{K}' \right) \mathbf{I}_{\rm sun} +
      (\Psi_M+\Psi_0) \mathbf{S} \right],
      \label{eq:slab_delo}

   where the coefficients :math:`\Psi_M` and :math:`\Psi_0` depend only
   on the optical thickness of the slab at the frequency and
   line-of-sight under consideration, since their expressions are:

   .. math::

      \begin{aligned}
      \Psi_M&=& \frac{1-e^{-\tau}}{\tau} - e^{-\tau},\nonumber \\
      \Psi_0 &=&1-\frac{1-e^{-\tau}}{\tau}.\end{aligned}

   Note that Eq. ([eq:slab\_delo]) for the emergent Stokes vector is the
   one used by :raw-latex:`\cite{trujillo_asensio07}` for investigating
   the impact of atomic level polarization on the Stokes profiles of the
   He i 10830 Å multiplet. We point out that, strictly speaking, it can
   be considered only as the exact analytical solution of the
   optically-thin constant-property slab problem [3]_. The reason why
   Eq. ([eq:slab\_delo]) is, in general, an approximate expression for
   calculating the emergent Stokes vector is because its derivation
   assumes that the Stokes vector within the slab varies linearly with
   the optical distance. However, it provides a fairly good
   approximation to the emergent Stokes profiles (at least for all the
   problems we have investigated in this paper). Moreover, the results
   of fig. 2 of :raw-latex:`\cite{trujillo_asensio07}` remain also
   virtually the same when using instead the exact Eq. ([eq:slab\_peo]),
   which from a computational viewpoint is significantly less efficient
   than the approximate Eq. ([eq:slab\_delo]).

-  *Exact analytical solution of the problem of a constant-property slab
   when neglecting the second-order terms of the Stokes-vector transfer
   equation*. Simplified expressions for the emergent Stokes vector can
   be obtained when :math:`\epsilon_I{\gg}\epsilon_X` and
   :math:`\eta_I{\gg}(\eta_X,\rho_X)`, which justifies to neglect the
   second-order terms of Eq. ([eq:rad\_transfer]). The resulting
   approximate formulae for the emergent Stokes parameters are given by
   eqs. (9) and (10) of :raw-latex:`\cite{trujillo_asensio07}`, which
   are identical to those used by :raw-latex:`\cite{trujillo_merenda05}`
   for modeling the Stokes profiles observed in solar chromospheric
   spicules. We point out that there is a typing error in the sentence
   that introduces such eqs. (9) and (10) in
   :raw-latex:`\cite{trujillo_asensio07}`, since they are obtained only
   when the above-mentioned second-order terms are neglected in Eq.
   ([eq:rad\_transfer]), although it is true that there are no
   magneto-optical terms in the resulting equations.

-  *Optically thin limit*. Finally, the most simple solution is obtained
   when taking the optically thin limit (:math:`\tau{\ll}1`) in the
   equations reported in the previous point, which lead to the equations
   (11) and (12) of :raw-latex:`\cite{trujillo_asensio07}`. Note that if
   :math:`\mathbf{I}_{\rm sun}=0` (i.e., :math:`I_0=X_0=0`), then such
   optically thin equations imply that
   :math:`{X/I}\,{\approx}\,{\epsilon_X}/{\epsilon_I}`.

The coefficients of the emission vector and of the propagation matrix
depend on the multipolar components, :math:`\rho^K_Q(J,J^{'})`, of the
atomic density matrix. Let us recall now the meaning of these physical
quantities and how to calculate them in the presence of an arbitrary
magnetic field under given illumination conditions.

The multipolar components of the atomic density matrix
------------------------------------------------------

We quantify the atomic polarization of the atomic levels using the
multipolar components of the atomic density matrix. We assume that the
atom can be correctly described under the framework of the
:math:`L`-:math:`S` coupling
:raw-latex:`\citep[e.g.,][]{condon_shortley35}`. The different
:math:`J`-levels are grouped in terms with well defined values of the
electronic angular momentum :math:`L` and the spin :math:`S`. We neglect
the influence of hyperfine structure and assume that the energy
separation between the :math:`J`-levels pertaining to each term is very
small in comparison with the energy difference between different terms.
Therefore, we allow for coherences between different :math:`J`-levels
pertaining to the same term but not between the :math:`J`-levels
pertaining to different terms. As a result, we can represent the atom
under the formalism of the multi-term atom discussed by
:raw-latex:`\cite{landi_landolfi04}`.

In the absence of magnetic fields the energy eigenvectors can be written
using Dirac’s notation as :math:`|\beta L S J M\rangle`, where
:math:`\beta` indicates a set of inner quantum numbers specifying the
electronic configuration. In general, if a magnetic field of arbitrary
strength is present, the vectors :math:`|\beta L S J M\rangle` are no
longer eigenfunctions of the total Hamiltonian and :math:`J` is no
longer a good quantum number. In this case, the eigenfunctions of the
full Hamiltonian can be written as the following linear combination:

.. math::

   \label{eq:eigenfunctions_total_hamiltonian}
   |\beta L S j M\rangle = \sum_J C_J^j(\beta L S, M) |\beta L S J M\rangle,

 where :math:`j` is a pseudo-quantum number which is used for labeling
the energy eigenstates belonging to the subspace corresponding to
assigned values of the quantum numbers :math:`\beta`, :math:`L`,
:math:`S`, and :math:`M`, and where the coefficients :math:`C_J^j` can
be chosen to be real.

In the presence of a magnetic field sufficiently weak so that the
magnetic energy is much smaller than the energy intervals between the
:math:`J`-levels, the energy eigenvectors are still of the form
:math:`|\beta L S J M\rangle`
(:math:`C_J^j(\beta L S, M) \approx \delta_{Jj}`), and the splitting of
the magnetic sublevels pertaining to each :math:`J`-level is linear with
the magnetic field strength. For stronger magnetic fields, we enter the
incomplete Paschen-Back effect regime in which the energy eigenvectors
are of the general form given by Eq.
([eq:eigenfunctions\_total\_hamiltonian]), and the splitting among the
various :math:`M`-sublevels is no longer linear with the magnetic
strength. If the magnetic field strength is further increased we
eventually reach the so-called complete Paschen-Back effect regime,
where the energy eigenvectors are of the form
:math:`|L S M_L M_S\rangle` and each :math:`L`-:math:`S` term splits
into a number of components, each of which corresponding to particular
values of (:math:`M_L+2M_S`).

Within the framework of the multi-term atom model the atomic
polarization of the energy levels is described with the aid of the
density matrix elements

.. math:: \rho^{\beta L S}(jM,j'M') = \langle \beta L S j M | \rho | \beta L S j' M'\rangle,

 where :math:`\rho` is the atomic density matrix operator. Using the
expression of the eigenfunctions of the total Hamiltonian given by Eq.
([eq:eigenfunctions\_total\_hamiltonian]), the density matrix elements
can be rewritten as:

.. math::

   \rho^{\beta L S}(jM,j'M') = \sum_{JJ'} C_J^j(\beta L S, M) C_{J'}^{j'}(\beta L
   S, M') \rho^{\beta L S}(JM,J'M'),

 where :math:`\rho^{\beta L S}(JM,J'M')` are the density matrix elements
on the basis of the eigenvectors :math:`| \beta L S J M\rangle`.

Following :raw-latex:`\cite{landi_landolfi04}`, it is helpful to use the
spherical statistical tensor representation, which is related to the
previous one by the following linear combination:

.. math::

   \begin{aligned}
   {^{\beta LS}\rho^K_Q(J,J')} &=& \sum_{jj'MM'} C_J^j(\beta L S, M)
   C_{J'}^{j'}(\beta L S, M') \nonumber \\
   &\times& (-1)^{J-M} \sqrt{2K+1} { \left(\begin{array}{ccc}
   J&J'&K\\
   M&-M'&-Q
   \end{array}\right) } 
   \rho^{\beta L S}(jM,j'M'),\end{aligned}

 where the 3-j symbol is defined as indicated by any suitable textbook
on Racah algebra.

Statistical equilibrium equations
---------------------------------

In order to obtain the :math:`{^{\beta LS}\rho^K_Q(J,J')}` elements we
have to solve the statistical equilibrium equations. These equations,
written in a reference system in which the quantization axis (:math:`Z`)
is directed along the magnetic field vector and neglecting the influence
of collisions, can be written as :raw-latex:`\citep{landi_landolfi04}`:

.. math::

   \begin{aligned}
   \frac{d}{dt} {^{\beta LS}\rho^K_Q(J,J')} &=& -2\pi \mathrm{i} \sum_{K' Q'}
   \sum_{J'' J'''} N_{\beta LS}(KQJJ',K'Q'J''J''') {^{\beta LS}\rho^{K'}_{Q'}(J'',J''')}
   \nonumber \\
   &+& \sum_{\beta_\ell L_\ell K_\ell Q_\ell J_\ell J_\ell'} {^{\beta_\ell L_\ell
   S}\rho^{K_\ell}_{Q_\ell}(J_\ell,J_\ell')} 
   \mathbb{T}_A(\beta L S K Q J J', \beta_\ell L_\ell S K_\ell Q_\ell J_\ell
   J_\ell') \nonumber \\
   &+& \sum_{\beta_u L_u K_u Q_u J_u J_u'} {^{\beta_u L_u
   S}\rho^{K_u}_{Q_u}(J_u,J_u')} 
   \Big[ \mathbb{T}_E(\beta L S K Q J J', \beta_u L_u S K_u Q_u J_u J_u') \nonumber \\
   & &\qquad \qquad \qquad \qquad \qquad + \mathbb{T}_S(\beta L S K Q
   J J', \beta_u L_u S K_u Q_u J_u J_u') \Big] \nonumber \\
   &-& \sum_{K' Q' J'' J'''} {^{\beta L S}\rho^{K'}_{Q'}(J'',J''') } \Big[
   \mathbb{R}_A(\beta L S K Q J J' K' Q' J'' J''') \nonumber \\
   & & + \mathbb{R}_E(\beta L S K Q J J' K'
   Q' J'' J''') + \mathbb{R}_S(\beta L S K Q J J' K' Q' J'' J''') \Big].
   \label{eq:see}\end{aligned}

 The first term in the right hand side of Eq. ([eq:see]) takes into
account the influence of the magnetic field on the atomic level
polarization. This term has its simplest expression in the chosen
magnetic field reference frame
:raw-latex:`\citep[see eq. 7.41 of][]{landi_landolfi04}`. In any other
reference system, a more complicated expression arises. The second,
third and fourth terms account, respectively, for coherence transfer due
to absorption from lower levels (:math:`\mathbb{T}_A`), spontaneous
emission from upper levels (:math:`\mathbb{T}_E`) and stimulated
emission from upper levels (:math:`\mathbb{T}_S`). The remaining terms
account for the relaxation of coherences due to absorption to upper
levels (:math:`\mathbb{R}_A`), spontaneous emission to lower levels
(:math:`\mathbb{R}_E`) and stimulated emission to lower levels
(:math:`\mathbb{R}_S`), respectively.

The stimulated emission and absorption transfer and relaxation rates
depend explicitly on the radiation field properties
:raw-latex:`\citep[see eqs. 7.45 and 7.46 of][]{landi_landolfi04}`. The
symmetry properties of the radiation field are accounted for by the
spherical components of the radiation field tensor:

.. math::

   J^K_Q(\nu) = \oint \frac{d\Omega}{4\pi} \sum_{i=0}^3
   \mathcal{T}^K_Q(i,\mathbf{\Omega}) S_i(\nu,\mathbf{\Omega}).
   \label{eq:jkq}

The quantities :math:`\mathcal{T}^K_Q(i,\mathbf{\Omega})` are spherical
tensors that depend on the reference frame and on the ray direction
:math:`\mathbf{\Omega}`. They are given by

.. math::

   \mathcal{T}^K_Q(i,\mathbf{\Omega}) = \sum_P t^K_P(i) \mathcal{D}^K_{PQ}(R'),
   \label{eq:tkq}

 where :math:`R'` is the rotation that carries the reference system
defined by the line-of-sight :math:`\mathbf{\Omega}` and by the
polarization unit vectors :math:`\mathbf{e}_1` and :math:`\mathbf{e}_2`
into the reference system of the magnetic field, while
:math:`\mathcal{D}^K_{PQ}(R')` is the usual rotation matrix
:raw-latex:`\citep[e.g.,][]{edmonds60}`. Table 5.6 in
:raw-latex:`\cite{landi_landolfi04}` gives the
:math:`\mathcal{T}^K_Q(i,\mathbf{\Omega})` values for each Stokes
parameter :math:`S_i` (with :math:`S_0=I`, :math:`S_1=Q`, :math:`S_2=U`
and :math:`S_3=V`).

Emission and absorption coefficients
------------------------------------

Once the multipolar components :math:`{^{\beta L S}\rho^{K}_{Q}(J,J') }`
are known, the coefficients :math:`\epsilon_I` and :math:`\epsilon_X`
(with :math:`X=Q,U,V`) of the emission vector and the coefficients
:math:`\eta_I`, :math:`\eta_X`, and :math:`\rho_X` of the propagation
matrix for a given transition between an upper term
:math:`(\beta L_u S)` and an lower term :math:`(\beta L_\ell S)` can be
calculated with the expressions of §7.6.b in
:raw-latex:`\cite{landi_landolfi04}`. These radiative transfer
coefficients are proportional to the number density of atoms,
:math:`\mathcal{N}`. Their defining expressions contain also the Voigt
profile and the Faraday-Voigt profile
:raw-latex:`\citep[see \S5.4 in][]{landi_landolfi04}`, which involve the
following parameters: :math:`a` (i.e., the reduced damping constant),
:math:`v_\mathrm{th}` (i.e., the velocity that characterizes the thermal
motions, which broaden the line profiles), and :math:`v_\mathrm{mac}`
(i.e., the velocity of possible bulk motions in the plasma, which
produce a Doppler shift).

It is important to emphasize that the expressions for the emission and
absorption coefficients and those of the statistical equilibrium
equations are written in the reference system whose quantization axis is
parallel to the magnetic field. The following equation indicates how to
obtain the density matrix elements in a new reference system:

.. math::

   \left[ {^{\beta L S}\rho^{K}_{Q}(J,J') } \right]_\mathrm{new} = \sum_{Q'} \left[
   {^{\beta L S}\rho^{K}_{Q'}(J,J') } \right]_\mathrm{old}
   \mathcal{D}^K_{Q' Q}(R)^*,

 where :math:`\mathcal{D}^K_{Q' Q}(R)^*` is the complex conjugate of the
rotation matrix for the rotation :math:`R` that carries the old
reference system into the new one.

Inversion
=========

Our inversion strategy is based on the minimization of a merit function
that quantifies how well the Stokes profiles calculated in our
atmospheric model reproduce the observed Stokes profiles. To this end,
we have chosen the standard :math:`\chi^2`–function, defined as:

.. math::

   \chi^2 = \frac{1}{4N_\lambda} \sum_{i=1}^4 \sum_{j=1}^{N_\lambda} 
   \frac{\left[S_i^\mathrm{syn}(\lambda_j)-S_i^\mathrm{obs}(\lambda_j) \right]^2}{
   \sigma_i^2(\lambda_j)} ,

 where :math:`N_\lambda` is the number of wavelength points and
:math:`\sigma_i^2(\lambda_j)` is the variance associated to the
:math:`j`-th wavelength point of the :math:`i`-th Stokes profiles. The
minimization algorithm tries to find the value of the parameters of our
model that lead to synthetic Stokes profiles :math:`S_i^\mathrm{syn}`
with the best possible fit to the observations. For our slab model, the
number of parameters (number of dimensions of the :math:`\chi^2`
hypersurface) lies between 5 and 7, the maximum value corresponding to
the optically thick case. The magnetic field vector (:math:`B`,
:math:`\theta_B` and :math:`\chi_B`), the thermal velocity
(:math:`v_\mathrm{th}`) and the macroscopic velocity
(:math:`v_\mathrm{mac}`) are always required. This set of parameters is
enough for the case of an optically thin slab. In order to account for
radiative transfer effects, we need to define the optical depth of the
slab along its normal direction and at a suitable reference wavelength
(e.g., the central wavelength of the red blended component for the 10830
Å multiplet). In addition, we may additionally need to include the
damping parameter (:math:`a`) of the Voigt profile if the wings of the
observed Stokes profiles cannot be fitted using Gaussian line profiles.

Global Optimization techniques
------------------------------

In order to avoid the possibility of getting trapped in a local minimum
of the :math:`\chi^2` hypersurface, global optimization methods have to
be used. We have chosen the DIRECT algorithm
:raw-latex:`\citep{Jones_DIRECT93}`, whose name derives from one of its
main features: *di*\ viding *rect*\ angles. The idea is to recursively
sample parts of the space of parameters, improving in each iteration the
location of the part of the space where the global minimum is
potentially located. The decision algorithm is based on the assumption
that the function is Lipschitz continuous
:raw-latex:`\citep[see][for details]{Jones_DIRECT93}`. The method works
very well in practice and can indeed find the minimum in functions that
do not fulfill the condition of Lipschitz continuity. The reason is that
the DIRECT algorithm does not require the explicit calculation of the
Lipschitz constant but it uses all possible values of such a constant to
determine if a region of the parameter space should be broken into
subregions because of its potential interest
:raw-latex:`\citep[see][for details]{Jones_DIRECT93}`.

Since the intensity profile is not very sensitive to the presence of a
magnetic field (at least for magnetic field strengths of the order of or
smaller than 1000 G), we have decided to estimate the optical thickness
of the slab, the thermal and the macroscopic velocity of the plasma and
the damping constant by using only the Stokes :math:`I` profile, and
then to determine the magnetic field vector by using the polarization
profiles. The full inversion scheme begins by applying the DIRECT method
to obtain a first estimation of the indicated four parameters by using
only Stokes :math:`I`. Afterwards, some LM iterations are carried out to
refine the initial values of the model’s parameters obtained in the
previous step. Once the LM method has converged, the inferred values of
:math:`v_\mathrm{th}`, :math:`v_\mathrm{mac}` (together with :math:`a`
and :math:`\Delta \tau`, when these are parameters of the model) are
kept fixed in the next steps, in which the DIRECT method is used again
for obtaining an initial approximation of the magnetic field vector
(:math:`B`,\ :math:`\theta_B`,\ :math:`\chi_B`). According to our
experience, the first estimate of the magnetic field vector given by the
DIRECT algorithm is typically very close to the final solution.
Nevertheless, some iterations of the LM method are performed to refine
the value of the magnetic field strength, inclination and azimuth. In
any case, although we have found very good results with this procedure,
the specific inversion scheme is fully configurable and can be tuned for
specific problems.

Our experience has proved that the following strategy is appropriate for
inverting prominences. Two initial DIRECT+LM cycles with weights
:math:`(1,0,0,0)` to invert the thermodynamical parameters. Then, two
DIRECT+LM cycles in which :math:`B`, :math:`\theta_B` and :math:`\chi_B`
are left free with weights :math:`(0,0.1,0.1,1)` which tries to set the
correct polarity of the field given by Stokes :math:`V`. An additional
LM cycle in which we fit only :math:`\theta_B` and :math:`\chi_B` with
the weights :math:`(0,1,1,0.3)` and a last LM cycle with weights
:math:`(0,0.3,0.3,1)` leaving the full magnetic field vector free.

Convergence
-----------

We let the DIRECT algorithm locate the global minimum in a region whose
hypervolume is :math:`V`. This hypervolume is obtained as the product of
the length :math:`d_i` of each dimension associated with each of the
:math:`N` parameters:

.. math:: V = \prod_i^N d_i.

 When the hypervolume decreases by a factor :math:`f` after the DIRECT
algorithm has discarded some of the hyperrectangles, its size along each
dimension is approximately decreased by a factor :math:`f^{1/N}`. In
order to end up with a small region where the global minimum is located,
many subdivisions are necessary, thus requiring many function
evaluations.

The most time consuming part of any optimization procedure is the
evaluation of the merit function. The DIRECT algorithm needs only a
reduced number of evaluations of the merit function to find the region
where the global minimum is located. For this reason, we have chosen it
as the initialization part of the LM method. Since the initialization
point is close to the global minimum, the LM method, thanks to its
quadratic behavior, rapidly converges to the minimum.

Stopping criterium
------------------

We have used two stopping criteria for the DIRECT algorithm. The first
one is stopping when the ratio between the hypervolume where the global
minimum is located and the original hypervolume is smaller than a given
threshold. This method has been chosen when using the DIRECT algorithm
as an initialization for the LM method, giving very good results. The
other good option, suggested by :raw-latex:`\cite{Jones_DIRECT93}`, is
to stop after a fixed number of evaluations of the merit function.

Ambiguities in the Hanle effect in the saturation regime
========================================================

In the saturation regime of the Hanle effect, Stokes :math:`Q` and
:math:`U` are insensitive to the field strength, but are sensitive to
the geometry of the field. For a :math:`J=0 \to J=1` transition, the
linear polarization can be written as:

.. math::

   \begin{aligned}
   Q &=& \frac{q}{2} \left( 3 \cos^2 \theta_B-1 \right) \sin^2\Theta_B \cos 2\Phi_B \nonumber \\
   U &=& \frac{q}{2} \left( 3 \cos^2 \theta_B-1 \right) \sin^2\Theta_B \sin 2\Phi_B.\end{aligned}

 These expressions contain a mixture of angles to make it clear that the
polarization amplitude depends on both the angle between the vertical
and the magnetic field and between the magnetic field and the
line-of-sight (LOS).

The coordinates of the magnetic field vector :math:`\mathbf{B}` in the
reference system of the vertical and the reference system of the LOS
are:

.. math::

   \begin{aligned}
   \mathbf{B} &=& B \left(\sin \theta_B \cos \phi_B \mathbf{i}+\sin \theta_B \sin \phi_B \mathbf{j}+\cos \theta_B \mathbf{k} \right) \nonumber \\
   \mathbf{B} &=& B \left(\sin \Theta_B \cos \Phi_B \mathbf{i}'+\sin \Theta_B \sin \Phi_B \mathbf{j}'+\cos \Theta_B \mathbf{k}' \right),\end{aligned}

 where the unit vectors are related by a simple rotation:

.. math::

   \begin{aligned}
   \mathbf{i}' &=& \cos \theta \mathbf{i} - \sin \theta \mathbf{k} \nonumber \\
   \mathbf{k}' &=& \sin \theta \mathbf{i} + \cos \theta \mathbf{k}.\end{aligned}

 Introducing these relations on the expression for the magnetic field,
we find that the following has to be fulfilled, given that the magnetic
field vector is the same in both reference systems:

.. math::

   \begin{aligned}
   \sin \theta_B \cos \phi_B &=& \sin \Theta_B \cos \Phi_B \cos \theta + \cos \Theta_B + \sin \theta \nonumber \\
   \sin \theta_B \sin \phi_B &=& \sin \Theta_B \sin \Phi_B \nonumber \\
   \cos \theta_B &=& \cos \Theta_B \cos \theta - \sin \Theta_B \cos \Phi_B \sin \theta.\end{aligned}

 Solving the previous three equations in the two directions, we find the
following transformations between the angles in the vertical reference
system and the LOS reference system:

.. math::

   \begin{aligned}
   \cos \Theta_B &=& \cos\theta \cos\theta_B + \sin\theta \sin\theta_B \cos\phi_B \nonumber \\
   \sin \Theta_B &=& +\sqrt{1-\cos^2\Theta_B} \nonumber \\
   \cos \Phi_B &=& \frac{\cos\theta \sin\theta_B \cos\phi_B - \cos\theta_B \sin\theta}{\sin \Theta_B} \nonumber \\
   \sin \Phi_B &=& \frac{\sin\theta_B \sin\phi_B}{\sin\Theta_B}\end{aligned}

 and

.. math::

   \begin{aligned}
   \cos \theta_B &=& \cos\theta \cos\Theta_B - \sin\theta \sin\Theta_B \cos\Phi_B \nonumber \\
   \sin \theta_B &=& +\sqrt{1-\cos^2\theta_B} \nonumber \\
   \cos \phi_B &=& \frac{\cos\theta \sin\Theta_B \cos\Phi_B + \cos\Theta_B \sin\theta}{\sin \theta_B} \nonumber \\
   \sin \phi_B &=& \frac{\sin\Theta_B \sin\Phi_B}{\sin\theta_B}.\end{aligned}

 Note that, since :math:`\Theta_B \in [0,\pi]`, we can safely use the
square root and take the positive value. In order to transform from one
reference system to the other, we can compute the inclination easily by
inverting the sinus or the cosinus. However, the situation is different
for the azimuth, because the range of variation is :math:`[-\pi,\pi]`.
Therefore, one has to compute the cosinus and the sinus separately and
the decide which is the correct quadrant fo the angle in terms of the
signs of both quantities.

Four possible kinds of ambiguities can exist for the Stokes :math:`Q`
and :math:`U` parameters. The idea is that :math:`\Phi_B` can be
modified and still obtain the same :math:`Q` and :math:`U` by properly
adjusting the value of :math:`\Theta_B`. It is clear that, given that
the term that can be used to compensate for the change in the azimuth on
the LOS reference system is the same for Stokes :math:`Q` and :math:`U`,
we can only compensate for changes in the sign. Therefore, we have the
following potential ambiguities:

.. math::

   \begin{aligned}
   \Phi_B' &=& \Phi_B \nonumber \\
   \Phi_B' &=& \Phi_B -\pi/2 \nonumber \\
   \Phi_B' &=& \Phi_B + \pi/2 \nonumber \\
   \Phi_B' &=& \Phi_B + \pi.\end{aligned}

 For each case, we have to compute the value of :math:`\Theta_B'` that
keeps the value of :math:`Q` and :math:`U` unchanged. Therefore, once we
find a solution to the inversion problem in the form of the pair
:math:`(\theta_B,\phi_B)`, we can find the remaining solutions in the
saturation regime following the recipes that we present now. Remember
that, unless one knows the polarity of the field, or in other words, the
sign :math:`\cos\Theta_B`, the number of potential ambiguous solutions
is 8. If the polarity of the field is known, the number is typically
reduced to 4 (or 2 if no 90\ :math:`^\circ` ambiguity is present).

PhiB’=PhiB
==========

Under this change, we have that

.. math:: \cos 2\Phi_B' = \cos 2\Phi_B, \quad \sin 2\Phi_B' = \sin 2\Phi_B, \quad \cos \Phi_B' = \cos \Phi_B, \quad \sin \Phi_B' = \sin \Phi_B.

 Making use of the previous relations between the angles wrt to the
vertical and the LOS, we have to solve the following equation:

.. math:: \left( 3 \cos^2\theta_B'-1 \right) \sin^2 \Theta_B' = \left( 3 \cos^2\theta_B-1 \right) \sin^2 \Theta_B,

 which can be written as:

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta - \sin\theta \sin\Theta_B' \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

 After some algebra and doing the substitution :math:`t=\sin\Theta_B'`,
we end up with the following equation to be solved:

.. math:: A t^4 + Bt^2 + C t^3 \sqrt{1-t^2} = K,

 where

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \cos^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& -6 \cos\theta \sin\theta \cos \Phi_B \nonumber \\
   K &=& \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

 The previous equation can be solved if we make the change of variables
:math:`t=\pm \sqrt{Z}`, resulting in:

.. math:: (C^2+A^2) Z^4 + (-C^2+2AB) Z^3 + (-2AK+B^2) Z^2 - 2BKZ + K^2 = 0.

 This polynomial of 4-th order can have four different solutions. From
these solutions, we have to take only the real solutions which are
larger than 0, given the range of variation of :math:`\Theta_B`:

.. math:: t \in \mathbb{R}, \qquad 0 \leq t \leq 1.

 Once the solutions for :math:`t` are found, we make
:math:`\Theta_B' = \arcsin t`. Note that, for a fixed value of
:math:`t`, two values of :math:`\Theta_B'` are possible. We choose the
correct one by evaluating the expressions for :math:`Q` and :math:`U`
and testing which of the two possible choices give the values equal (or
very similar) to the original ones.

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B)` to the vertical reference
system.

PhiB’=PhiB+pi
=============

Under this change, we have:

.. math:: \cos 2\Phi_B' = \cos 2\Phi_B, \quad \sin 2\Phi_B' = \sin 2\Phi_B, \quad \cos \Phi_B' = -\cos \Phi_B, \quad \sin \Phi_B' = -\sin \Phi_B.

 Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

 The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \cos^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& 6 \cos\theta \sin\theta \cos \Phi_B \nonumber \\
   K &=& \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B+\pi)` to the vertical
reference system.

PhiB’=PhiB+pi/2
===============

Under this change, we have:

.. math:: \cos 2\Phi_B' = -\cos 2\Phi_B, \quad \sin 2\Phi_B' = -\sin 2\Phi_B, \quad \cos \Phi_B' = -\sin \Phi_B, \quad \sin \Phi_B' = \cos \Phi_B.

 Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \sin\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

 The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \sin^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& 6 \cos\theta \sin\theta \sin \Phi_B \nonumber \\
   K &=& -\left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B+\pi/2)` to the vertical
reference system.

PhiB’=PhiB-pi/2
===============

Under this change, we have:

.. math:: \cos 2\Phi_B' = -\cos 2\Phi_B, \quad \sin 2\Phi_B' = -\sin 2\Phi_B, \quad \cos \Phi_B' = \sin \Phi_B, \quad \sin \Phi_B' = -\cos \Phi_B.

 Following the same approach, we have to solve for :math:`\Theta_B'` in

.. math::

   \left[ 3 \left( \cos \Theta_B' \cos \theta + \sin\theta \sin\Theta_B' \sin\Phi_B\right)^2-1 \right] \sin^2 \Theta_B' = 
   \left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.

 The solution are obtained as the roots of the same equations as before
but now

.. math::

   \begin{aligned}
   A &=& -3\cos^2 \theta + 3\sin^2 \theta \sin^2 \Phi_B \nonumber \\
   B &=& 3\cos^2 \theta - 1 \nonumber \\
   C &=& -6 \cos\theta \sin\theta \sin \Phi_B \nonumber \\
   K &=& -\left[ 3 \left( \cos \Theta_B \cos \theta - \sin\theta \sin\Theta_B \cos\Phi_B\right)^2-1 \right] \sin^2 \Theta_B.\end{aligned}

The angles :math:`(\theta_B,\phi_B)` are obtained by doing the
transformation from :math:`(\Theta_B',\Phi_B-\pi/2)` to the vertical
reference system.


11 natexlab#1#1

.. [Ref]_, L., Trujillo Bueno, J., & Landi Degl’Innocenti, E. 2007, ApJ, 666, 588

, E. U., & Shortley, G. H. 1935, The Theory of Atomic Spectra
(Cambridge: Cambridge University Press)

, A. R. 1960, Angular Momentum in Quantum Mechanics (Princeton
University Press)

, D. R., Perttunen, C. D., & Stuckmann, B. E. 1993, Journal of
Optimization Theory and Applications, 79, 157

, E., & Landi Deglinnocenti, M. 1985, Sol. Phys., 97, 239

, E., & Landolfi, M. 2004, Polarization in Spectral Lines (Kluwer
Academic Publishers)

, K. 2000, in Allen’s Astrophysical Quantities, ed. A. N. Cox (New York:
Springer Verlag and AIP Press)

, D. E., Durrant, C. J., & Murphy, G. A. 1989, ApJ, 339, 1093

, J. 2003, in Stellar Atmosphere Modeling, ed. I. Hubeny, D. Mihalas, &
K. Werner, ASP Conf. Ser. 288 (San Francisco: ASP), 551

, J., & Asensio Ramos, A. 2007, ApJ, 655, 642

, J., Merenda, L., Centeno, R., Collados, M., & Landi Degl’Innocenti, E.
2005, ApJ, 619, L191

.. [1]
    (an acronym for HAnle and ZEeman Light) is one of the IAC computer
   programs for the synthesis and inversion of Stokes profiles resulting
   from the joint action of the Hanle and Zeeman effects.

.. [2]
   ``http://www.unidata.ucar.edu/software/netcdf/``

.. [3]
   More precisely, when the optical thickness of the slab is small in
   comparison with the eigenvalues of the matrix :math:`\mathbf{K}'`.

.. |Screen dump of the graphical front-end used for the inversion. [fig:inversion\_GUI]| image:: inv2.eps
.. |Screen dump of the graphical front-end used for the inversion. [fig:inversion\_GUI]| image:: inv3.eps

