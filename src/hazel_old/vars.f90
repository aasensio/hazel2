module vars
implicit none

	real(kind=8), parameter :: PK = 1.3806503d-16, UMA = 1.66053873d-24, PC = 2.99792458d10
   real(kind=8), parameter :: PH = 6.62606876d-27, PHK = PH / PK, PHC = 2.d0 * PH / PC**2
   real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10, PI = 3.14159265359d0
   real(kind=8), parameter :: PHC2 = 2.d0 * PH * PC**2, OPA = PI * PE**2 / (PME * PC)
   real(kind=8), parameter :: SQRTPI = 1.77245385091d0, EARTH_SUN_DISTANCE = 1.495979d13
   real(kind=8), parameter :: LARMOR = PE/(4.d0*PI*PME*PC), PMP = 1.67262158d-24, RSUN = 976.6d0
	
	
	integer :: isti, idep, imag, linear_solver, use_paschen_back, verbose_mode, working_mode, synthesis_mode
	integer :: use_mag_opt_RT, use_stim_emission_RT
	real(kind=8) :: delta_collision
	
	real(kind=8) :: fact(0:301)
	character(len=120) :: input_model_file, input_experiment, output_rho_vertical_upper, output_rho_magnetic_upper
	character(len=120) :: output_rho_vertical_lower, output_rho_magnetic_lower, output_rtcoef, output_rtcoef_zeeman
	character(len=120) :: input_observed_profiles, input_inverted_parameters, output_inverted_profiles, output_error_parameters
	character(len=120) :: direct_ranges, output_final_parameters
	integer :: n_terms, nrhos, file_pointer, ntran, is2, jlimit2, nparam
	
	real(kind=8), allocatable :: lsto2(:), energy(:,:)
	real(kind=8), allocatable :: epsilon(:,:), eta(:,:), epsilon_zeeman(:,:), eta_zeeman(:,:)
	real(kind=8), allocatable :: eta_stim(:,:), eta_stim_zeeman(:,:), mag_opt(:,:), mag_opt_zeeman(:,:)
	real(kind=8), allocatable :: mag_opt_stim(:,:), mag_opt_stim_zeeman(:,:)
	
	real(kind=8), allocatable :: threej(:)
	
!	real(kind=8) :: bgauss, thetabd, chibd, thb, chb, thetad, chid, gammad
!	real(kind=8) :: omin, omax, wl, vdopp
!	integer :: nemiss, no
	
	real(kind=8) :: dr(-1:1,-1:1), di(-1:1,-1:1)
	
	real(kind=8), allocatable :: SEE_A(:,:), SEE_b(:), SEE_mag_A(:,:)
	integer, allocatable :: ntab(:), ktab(:), qtab(:), irtab(:), j2tab(:), jp2tab(:), ntlsto(:), ntusto(:)
	real(kind=8), allocatable :: rnutab(:), aesto(:)
	
	real(kind=8), allocatable :: allen_ic(:,:), allen_cl(:,:)
	
	complex(kind=8), allocatable :: rhol(:,:,:,:), rhou(:,:,:,:), rhoml(:,:,:,:), rhomu(:,:,:,:)
	
	type variable_parameters
		real(kind=8) :: bgauss, thetabd, chibd, vdopp, dtau, delta_collision, vmacro, damping, beta, height, vdopp2, beta2
		real(kind=8) :: dtau2, vmacro2, bgauss2, thetabd2, chibd2, ff
		integer :: n_inverted, n_total, nslabs
		integer, pointer :: inverted(:)
	end type variable_parameters
	
	type fixed_parameters
		real(kind=8) :: thetad, chid, gammad, omin, omax, wl
		real(kind=8) :: Stokes_incident(0:3), Aul, Bul, Blu, nu
		integer :: no, nemiss, use_atomic_pol, total_forward_modeling
		integer :: pix_syn_id, col_syn_id, nlambda_syn_id, lambda_syn_id, map_syn_id, syn_id
		integer :: pix_par_id, col_par_id, map_par_id, par_id
		integer :: pix_error_id, col_error_id, map_error_id, error_id
		integer :: damping_treatment
		real(kind=8), pointer :: upper_direct(:), lower_direct(:), stokes_boundary(:,:)
		real(kind=8) :: volper
		integer :: DIRmaxf, stokes_boundary_len
		character(len=10) :: Stokes_incident_mode
		character(len=100) :: Stokes_incident_file
	end type fixed_parameters
	
	type type_observation
		integer :: n, npixel, ny, pix_id, col_id, nlambda_id, lambda_id, map_id, obs_id, nstokespar_id
		character(len=4) :: normalization
		integer :: observation_format, boundary_id, height_id, obstheta_id, obsgamma_id, parsInit_id, normalization_id
		real(kind=8), pointer :: wl(:), stokes(:,:), sigma(:,:)
	end type type_observation
	
	type type_inversion
		integer :: nf, iter_max, n_cycles, loop_cycle
		real(kind=8) :: lambda, chisq, chisq_old, min_chisq
		integer, pointer :: cycles(:,:), algorithm(:)
		real(kind=8), pointer :: stokes_weights(:,:)
		real(kind=8), pointer :: stokes_unperturbed(:,:), stokes_perturbed(:,:), dydx(:,:,:)
	end type type_inversion

	type atom_model
		integer :: ntran
		integer, pointer :: nterml(:), ntermu(:)
		real(kind=8), pointer :: ae(:), wavelength(:), reduction_factor(:), reduction_factor_omega(:), j10(:)
	end type atom_model
		
	
	type(variable_parameters) :: params, trial, scaled_params, errorparams
	type(fixed_parameters) :: fixed
	type(type_observation) :: observation
	type(type_inversion) :: inversion
	type(atom_model) :: atom
	
	character(len=6) :: parameters_name(18) = (/ 'B     ','thetaB','chiB  ', 'vdopp ', &
		'dtau  ','D^(2) ','v_mac ','damp  ','beta  ','h     ','dtau2 ', 'v_mac2',&
		'B2    ', 'thetB2', 'chiB2 ', 'vdopp2','ff1   ','beta2 '/)
	real(kind=8), parameter :: minim_pikaia(10) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -15.d0, 0.d0, 0.d0, 0.d0/)
	real(kind=8), parameter :: maxim_pikaia(10) = (/4000.d0, 180.d0, 180.d0, 20.d0, 3.d0, 18.d0, 40.d0, 10.d0, 10.d0, 100.d0/)

	integer :: starting_pixel, final_pixel
	
	real(kind=8) :: chi2Level

	real(kind=8), dimension(4) :: nbarExternal, omegaExternal
		
	
end module vars
