module pyHazelMod
use iso_c_binding, only: c_int, c_double
use vars
use maths
use io
use SEE
use rt_coef
use synth
use allen
implicit none

contains
subroutine c_hazel(index, B1Input, hInput, tau1Input, boundaryInput, &
    transInput, anglesInput, nLambdaInput, lambdaAxisInput, dopplerWidthInput, dampingInput, &
    dopplerVelocityInput, betaInput, nbarInput, omegaInput, &
    wavelengthOutput, stokesOutput) bind(c)

    integer(c_int), intent(in) :: transInput, index
    integer(c_int), intent(in) :: nLambdaInput
    real(c_double), intent(in), dimension(nLambdaInput) :: lambdaAxisInput
    real(c_double), intent(in), dimension(3) :: B1Input, anglesInput
    real(c_double), intent(in), dimension(4,nLambdaInput) :: boundaryInput
    real(c_double), intent(in), dimension(4) :: nbarInput, omegaInput
    real(c_double), intent(in) :: hInput, tau1Input, dopplerWidthInput, dampingInput, dopplerVelocityInput, betaInput
    real(c_double), intent(out), dimension(nLambdaInput) :: wavelengthOutput
    real(c_double), intent(out), dimension(4,nLambdaInput) :: stokesOutput

    integer :: n, nterml, ntermu, error
    
    real(c_double) :: ae, wavelength, reduction_factor, reduction_factor_omega, j10
    integer :: i, j
    logical :: recompute_see_rtcoef

    recompute_see_rtcoef = .True.

    if (dopplerVelocityInput_old == dopplerVelocityInput .and. dopplerWidthInput_old == dopplerWidthInput .and. dampingInput_old == dampingInput .and. all(B1Input_old == B1Input)) then
        recompute_see_rtcoef = .False.
    endif

    input_model_file = 'helium.mod'
    input_experiment = 'init_parameters.dat'        
    verbose_mode = 0
    linear_solver = 0       
    synthesis_mode = 5
    working_mode = 0

! Read the atomic model 
!   call read_model_file(input_model_file)
    
! Set the variables for the experiment from the parameters of the subroutine
    isti = 1
    imag = 1
    idep = 0
    use_paschen_back = 1
    params(index)%nslabs = 1
    
    params(index)%bgauss = B1Input(1)
    params(index)%thetabd = B1Input(2)
    params(index)%chibd = B1Input(3)
    params(index)%height = hInput
    params(index)%dtau = tau1Input
    params(index)%beta = betaInput

    fixed(index)%nemiss = transInput
    fixed(index)%use_atomic_pol = 1
    fixed(index)%thetad = anglesInput(1)
    fixed(index)%chid = anglesInput(2)
    fixed(index)%gammad = anglesInput(3)
            
    params(index)%vdopp = dopplerWidthInput
    params(index)%damping = dampingInput
    fixed(index)%damping_treatment = 0
    
! Read the wavelength of the transition that we want to synthesize
    fixed(index)%wl = atom%wavelength(transInput)
    
! Set the values of nbar and omega in case they are given
    fixed(index)%nbarExternal = nbarInput
    fixed(index)%omegaExternal = omegaInput

! If nbar=0 or omega=0, use the numbers from Allen. If not, treat them as reduction factors
    do i = 1, atom%ntran
        if (nbarInput(i) == 0) then
            fixed(index)%nbarExternal(i) = 1.0
        endif
        if (omegaInput(i) == 0) then
            fixed(index)%omegaExternal(i) = 1.0
        endif
    enddo
            
    params%vmacro = dopplerVelocityInput
    
    use_mag_opt_RT = 1
    use_stim_emission_RT = 1

!*********************************
!** SYNTHESIS MODE
!*********************************  
    fixed(index)%no = nLambdaInput
    observation(index)%n = fixed(index)%no
    
    if (.not.associated(observation(index)%wl)) allocate(observation(index)%wl(observation(index)%n))
    if (.not.associated(inversion(index)%stokes_unperturbed)) allocate(inversion(index)%stokes_unperturbed(0:3,fixed(index)%no))
    if (.not.associated(fixed(index)%stokes_boundary)) allocate(fixed(index)%stokes_boundary(0:3,observation(index)%n))
    
    fixed(index)%stokes_boundary(0:3,:) = boundaryInput

    observation(index)%wl = lambdaAxisInput

    fixed(index)%omax = minval(lambdaAxisInput)
    fixed(index)%omin = maxval(lambdaAxisInput)

    if (recompute_see_rtcoef) then
        
        ! Fill the statistical equilibrium equations
        call fill_SEE(params(index), fixed(index), 1, error)

        ! If the solution of the SEE gives an error, return
        if (error == 1) return
                
        ! Calculate the absorption/emission coefficients for a given transition
        call calc_rt_coef(params(index), fixed(index), observation(index), 1)
    
    endif

! Do the synthesis
    call do_synthesis(params(index), fixed(index), observation(index), inversion(index)%stokes_unperturbed, error)
    
    do i = 1, 4
        stokesOutput(i,:) = inversion(index)%stokes_unperturbed(i-1,:)
    enddo
    
    wavelengthOutput = observation(index)%wl + fixed(index)%wl

    dopplerVelocityInput_old = dopplerVelocityInput
    dopplerWidthInput_old = dopplerWidthInput
    dampingInput_old = dampingInput
    B1Input_old = B1Input


    ! if (allocated(epsilon)) deallocate(epsilon)
	! if (allocated(epsilon_zeeman)) deallocate(epsilon_zeeman)
    ! if (allocated(eta)) deallocate(eta)
    ! if (allocated(eta_zeeman)) deallocate(eta_zeeman)
    ! if (allocated(eta_stim)) deallocate(eta_stim)
    ! if (allocated(eta_stim_zeeman)) deallocate(eta_stim_zeeman)
    ! if (allocated(mag_opt)) deallocate(mag_opt)
    ! if (allocated(mag_opt_zeeman)) deallocate(mag_opt_zeeman)
    ! if (allocated(mag_opt_stim)) deallocate(mag_opt_stim)
    ! if (allocated(mag_opt_stim_zeeman)) deallocate(mag_opt_stim_zeeman)
    
!   open(unit=31,file=input_model_file,action='write',status='replace')
!   close(31,status='delete')
        
end subroutine c_hazel

subroutine c_init() bind(c)
        
! Initialize the random number generator
    call random_seed

! Initialize Allen's data
    call read_allen_data
        
! Fill the factorial array
    call factrl
                
    input_model_file = 'helium.mod'
    
! Read the atomic model 
    call read_model_file(input_model_file)
    
end subroutine c_init

end module pyHazelMod