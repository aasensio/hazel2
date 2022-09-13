module pyHazelMod
use iso_c_binding, only: c_int, c_double,c_char !EDGAR: remember altering this when touching .pyx
use vars
use maths
use io
use SEE
use rt_coef
use synth
use allen
implicit none

contains
subroutine c_hazel(index, synMethInput,B1Input, hInput, tau1Input, boundaryInput, &
    transInput, anglesInput, nLambdaInput, lambdaAxisInput, dopplerWidthInput, dampingInput, &
    j10Input,dopplerVelocityInput, betaInput, nbarInput, omegaInput, &
    atompolInput,magoptInput,stimemInput,nocohInput,dcolInput, & !EDGAR: too verbose for just 4 numbers
    wavelengthOut, stokesOut, epsOut, etaOut, stimOut,error) bind(c)

    !EDGAR: no need to read ntransInput because ntrans already initialized in atom%ntran
    integer(c_int), intent(in) :: synMethInput,transInput, index !EDGAR: index runs from 1 to n_chromo
    integer(c_int), intent(in) :: nLambdaInput,atompolInput,magoptInput,stimemInput,nocohInput
    real(c_double), intent(in), dimension(nLambdaInput) :: lambdaAxisInput
    real(c_double), intent(in), dimension(3) :: B1Input, anglesInput
    real(c_double), intent(in), dimension(4,nLambdaInput) :: boundaryInput
    !EDGAR: dim shouuld be atom%ntran, not 4 (4 is only for Helium)
    real(c_double), intent(in), dimension(atom%ntran) :: nbarInput, omegaInput 
    real(c_double), intent(in), dimension(3)  :: dcolInput
    !EDGAR: add j10Input . in or inout to avoid losing memory between consecutive python calls??
    real(c_double), intent(in), dimension(atom%ntran) :: j10Input  
    real(c_double), intent(in) :: hInput, tau1Input, dopplerWidthInput, dampingInput, dopplerVelocityInput, betaInput 
    real(c_double), intent(out), dimension(nLambdaInput) :: wavelengthOut
    real(c_double), intent(out), dimension(4,nLambdaInput) :: stokesOut,epsOut!EDGAR: containers for opt coeffs
    real(c_double), intent(out), dimension(7,nLambdaInput) :: etaOut,stimOut !EDGAR: containers for opt coeffs
    integer(c_int), intent(out) :: error

    integer :: n, nterml, ntermu
    
    real(c_double) :: ae, wavelength, reduction_factor, reduction_factor_omega !, j10 !EDGAR remove j10
    integer :: i, j
    logical :: recompute_see_rtcoef

    error = 0
    error_code = 0
    

    params(index)%recompute_see_rtcoef = .True.
    ! If the parameters on which the RT coefficients depend on change, then recompute the coefficients
    if (params(index)%dopplerVelocityInput_old == dopplerVelocityInput .and. &
        params(index)%dopplerWidthInput_old == dopplerWidthInput .and. &
        params(index)%dampingInput_old == dampingInput .and. &
        fixed(index)%thetad_old == anglesInput(1) .and. &
        fixed(index)%chid_old == anglesInput(2) .and. &
        fixed(index)%gammad_old == anglesInput(3) .and. &
        all(params(index)%B1Input_old == B1Input)) then
            params(index)%recompute_see_rtcoef = .False.
    endif

    !input_model_file = 'helium.mod'
    !input_experiment = 'init_parameters.dat'        
    linear_solver = 0       
    synthesis_method = synMethInput
    working_mode = 0 !0: synthesis. 1: inversion
    ! Read the atomic model 
    ! call read_model_file(input_model_file)

    use_paschen_back = 1 
    
    !--------------------------------
    isti = 1  !stimulated emission term in SEE 
    imag = 1  !Hanle term in SEE 
    idep = 0  !depolarizing collision term in SEE

    if (dcolInput(1) /= 0.0 .or. dcolInput(2) /= 0.0 .or. dcolInput(3) /= 0.0) then
        idep = 1  
        params%delta_collision = dcolInput(1)   
        params%delta_collk1 = dcolInput(2)
        params%delta_collk2 = dcolInput(3)
    endif

    ! atompol cannot be -1 anymore to kill only anisotropy (only can be 0 or 1).
    !We can kill anisotropy for K=2 by setting j20f=0 but alignment in general does not dissapear 
    !only killing that anisotropy because orientation can be transformed into alignment in the SEE.
    !so, or we kill both alignment and orientation setting atompol=0 or we kill both anisotropies
    fixed(index)%use_atomic_pol = atompolInput !0(no atompol),1 (all atompol)    
    use_mag_opt_RT = magoptInput
    use_stim_emission_RT = stimemInput
 
    params(index)%nocoh = -99 !this absurd value makes level coherences to be fully added in SEE when nocohInput=0
    if (nocohInput /= 0) params(index)%nocoh = nocohInput  !gives the atom level in which cohs will be deactivated
   !--------------------------------

    params(index)%nslabs = 1  !EDGAR:this is not being used anymore    
    params(index)%bgauss = B1Input(1)
    params(index)%thetabd = B1Input(2)
    params(index)%chibd = B1Input(3)
    params(index)%height = hInput
    params(index)%dtau = tau1Input
    params(index)%beta = betaInput

    fixed(index)%nemiss = transInput
    fixed(index)%thetad = anglesInput(1)
    fixed(index)%chid = anglesInput(2)
    fixed(index)%gammad = anglesInput(3)
            
    params(index)%vdopp = dopplerWidthInput
    params(index)%damping = dampingInput
    fixed(index)%damping_treatment = 0

! Read the wavelength of the transition that we want to synthesize
    fixed(index)%wl = atom%wavelength(transInput)
    
! Set the values of nbar and omega in case they are given
    fixed(index)%nbarExternal = nbarInput  !set to ones (Allen nbars) in chromosphere.py
    fixed(index)%omegaExternal = omegaInput !set to zeroes (no anisotropy) in chromosphere.py

! EDGAR: If nbar=0 or omega=0, use the numbers from Allen. If not, treat them as reduction factors
! CAUTION: the following loop superseed previous initialization of allen parameters setting them to 1.0.
! We change now this, leaving the possibility of setting them from main python program or file.

    ! do i = 1, atom%ntran
    !     if (nbarInput(i) == 0) then
    !         fixed(index)%nbarExternal(i) = 1.0
    !     endif
    !     if (omegaInput(i) == 0) then
    !         fixed(index)%omegaExternal(i) = 1.0
    !     endif
    !     !CAUTION: the following two numbers are always read from atom model and must be 1.0 to avoid reductions
    !     !They are unnecessary and could be eliminated because nbarexternal and omegaexternal are 
    !     !already used as reduction factors if we change their value to be different than 1 or 0.
    !     !we are commenting their set up to 1.0 here becausse they are already entering as 1.0 from
    !     !the atom file.
    !     !atom%reduction_factor(i)=1.0  
    !     !atom%reduction_factor_omega(i)=1.0

    ! enddo


    params%vmacro = dopplerVelocityInput
    
    !EDGAR: we read atom file already with the init() routine before 
    !calling the actual routine so atom%j10 is already initialized.
    do i = 1, atom%ntran
        atom%j10(i)=j10Input(i)   !overwrite thevalue of the file with the python  
        if (verbose_mode > 0)print*,'j10 Python input:',atom%j10(i)
        !transOutput(i)=atom%wavelength(i)  !get the central wavelegnths and take them out to python
    enddo

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

    if (params(index)%recompute_see_rtcoef) then
        
        !EDGAR: in myhazel this block was inside do_synthesis and the opt coeffs 
        !had extra dimension for the slab number. They now belong to in_fixed/fixed var.

        ! Fill the statistical equilibrium equations
        call fill_SEE(params(index), fixed(index))

        ! If the solution of the SEE gives an error, return        
        if (error_code == 1) then
            error = 1
            return        
        endif
                
        ! Calculate the absorption/emission coefficients for a given transition
        !EDGAR: fixed%eps, fixed%eta, etc are now opt coeffs belonging to in_fixed parameters 
        !structure type in vars. To extract them out to python main, we have to collect them 
        !for each index into an array. Their way out fortran must be the present interface 
        !subroutine, hence pyx file must be modified and Mod.spectrum should also contain 
        !the final vector-like etas,epsilons,etc, such that one can retrieve them.
        !velocities are considered in rt_coef.f90 with parameter va.  
        call calc_rt_coef(params(index), fixed(index), observation(index))

        ! If the calculation of the RT coefficients gives and error, return
        if (error_code == 1) then
            error = 1
            return        
        endif
    
        !the output containers for the vector optical coeffs are defined here
        !inversion(index)%stokes_unperturbed(i_stokes,nu)  --> stokesOutput(i_stokes,nu)
        !inversion(index)%stokes_unperturbed(i_stokes,nu)  --> Output

    endif

! Do the synthesis
    call do_synthesis(params(index), fixed(index), observation(index), inversion(index)%stokes_unperturbed, error)

    ! If the synthesis gives an error, return
    if (error_code == 1) then
        error = 1
        return        
    endif
    
    do i = 1, 4
        stokesOut(i,:) = inversion(index)%stokes_unperturbed(i-1,:)
    enddo
    
    wavelengthOut = observation(index)%wl + fixed(index)%wl
!-----------------------------------------------------------------------------

    !EDGAR:OLD output in myhazel1.0. 
    !With this output we waste space and we loose info about the stimulated coefficients
    !   etaOutput(1,1,:) = eta(0,:) - eta_stim(0,:)
    !   etaOutput(2,2,:) = eta(0,:) - eta_stim(0,:)
    !   etaOutput(3,3,:) = eta(0,:) - eta_stim(0,:)
    !   etaOutput(4,4,:) = eta(0,:) - eta_stim(0,:)
    !   etaOutput(1,2,:) = eta(1,:) - eta_stim(1,:)
    !   etaOutput(2,1,:) = eta(1,:) - eta_stim(1,:)     
    !   etaOutput(1,3,:) = eta(2,:) - eta_stim(2,:)
    !   etaOutput(3,1,:) = eta(2,:) - eta_stim(2,:)
    !   etaOutput(1,4,:) = eta(3,:) - eta_stim(3,:)
    !   etaOutput(4,1,:) = eta(3,:) - eta_stim(3,:)
    !   etaOutput(2,3,:) = mag_opt(3,:) - mag_opt_stim(3,:)  
    !   etaOutput(3,2,:) = -(mag_opt(3,:) - mag_opt_stim(3,:)) !EDGAR: the signs were wrong!!!!
    !   etaOutput(2,4,:) = -(mag_opt(2,:) - mag_opt_stim(2,:))  !EDGAR: 1,2,3: Q,U,V
    !   etaOutput(4,2,:) = mag_opt(2,:) - mag_opt_stim(2,:)
    !   etaOutput(3,4,:) = mag_opt(1,:) - mag_opt_stim(1,:)
    !   etaOutput(4,3,:) = -(mag_opt(1,:) - mag_opt_stim(1,:))

    !The folowing output is more useful, more compact 
    !Requires set up of the atomicPolInput parameter.
    !Note that inner optical coeffs start in index 0 while Out optical coeffs in 1.
    if (fixed(index)%use_atomic_pol==0)then !when no atompol case, we extract only zeeman coefs.
        do i = 1, 4
            epsOut(i,:) = fixed(index)%epsilon_zeeman(i-1,:) 
            etaOut(i,:) = fixed(index)%eta_zeeman(i-1,:) !eta_I,eta_Q,eta_U,eta_V
            stimOut(i,:)=fixed(index)%eta_stim_zeeman(i-1,:)*use_stim_emission_RT !
        enddo
        do i = 5, 7
            etaOut(i,:) = fixed(index)%mag_opt_zeeman(i-4,:)  !store rho1(Q),rho2(U),rho3(V) in etaOut(5),(6) y (7)
            stimOut(i,:) = fixed(index)%mag_opt_stim_zeeman(i-4,:)*use_stim_emission_RT  !ro_QUV_stim
        enddo
    else
        do i = 1, 4
            epsOut(i,:) = fixed(index)%epsilon(i-1,:)
            etaOut(i,:) = fixed(index)%eta(i-1,:) !eta_I,eta_Q,eta_U,eta_V
            stimOut(i,:)=fixed(index)%eta_stim(i-1,:)*use_stim_emission_RT !
        enddo
        do i = 5, 7
            etaOut(i,:) = fixed(index)%mag_opt(i-4,:)  !store rho1(Q),rho2(U),rho3(V) in etaOut(5),(6) y (7)
            stimOut(i,:) = fixed(index)%mag_opt_stim(i-4,:)*use_stim_emission_RT  !ro_QUV_stim
        enddo
    
    endif

    !eta_i=eta^A_i - eta^S_i and idem for rho (rho_i=rho^A_i - rho^S_i)
    !eta_i(1:4)=etaOut(1:4) - stimOut(1:4)
    !rho_i(1:3)=etaOut(5:7) - stimOut(5:7)  
!-----------------------------------------------------------------------------

    params(index)%dopplerVelocityInput_old = dopplerVelocityInput
    params(index)%dopplerWidthInput_old = dopplerWidthInput
    params(index)%dampingInput_old = dampingInput
    params(index)%B1Input_old = B1Input
    fixed(index)%thetad_old = anglesInput(1)
    fixed(index)%chid_old = anglesInput(2)
    fixed(index)%gammad_old = anglesInput(3)


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
    RETURN
          
end subroutine c_hazel

PURE FUNCTION array_to_string(a)  RESULT (s)    ! copy char array to string
    !This func is needed to pass from array of chars(C way of dealing with strings) to fortran strings
    !Incompatibility issue will appear otherwise because we're interfacing using the ISO_C_BINDING 
    !features. The C language only defines strings as arrays of characters.  
    !Therefore, if receiving a string from C, it will be an array. 
    !To use this array as a string in Fortran, it has to be converted using a custom function like this
    CHARACTER,INTENT(IN) :: a(:)
    CHARACTER(SIZE(a)) :: s
    INTEGER :: i
    DO i = 1,SIZE(a)
       s(i:i) = a(i)
    END DO
END FUNCTION array_to_string


subroutine c_init(nchar,atomfileInput,verbose,ntransOut) bind(c)
    integer(c_int), intent(in) :: nchar  !EDGAR: passes the length of the next string
    character(c_char), intent(in) :: atomfileInput(nchar) !EDGAR:a C array of characters is entering
    integer(c_int), intent(in) ::verbose
    integer(c_int), intent(out) ::ntransOut
    integer :: i
    
    ntransOut=0

! Initialize the random number generator
    call random_seed

! Initialize Allen's data
    call read_allen_data
        
! Fill the factorial array
    call factrl    

    !NOW we use .atom files with relative paths
    !input_model_file = '../hazel/data/helium.atom'
    !input_model_file = '../hazel/data/sodium_hfs.atom'
    !convert from C array of chars to F90 string:
    input_model_file = '../hazel/data/'//array_to_string(atomfileInput)

    !EDGAR: #Not possible to call to logger.info from here.
    if (verbose > 0)print*,'Reading atom file:',input_model_file
    !EDGAR: here atom%j10(x) for each transition x is initialized to 0
    !and atom%ntran is read from file
    !c_init shouuld return atom%ntran to init nbar,j10 containers in python??
    call read_model_file_ok(input_model_file) 
    
    ntransOut=atom%ntran
    if (ntransOut==0)then
        print*,'No transitions in atom file. STOPPING.'
        STOP
    endif

    verbose_mode=verbose !EDGAR: verbose_mode is general fortran var, like input_model_file 

    ! Force recomputation of RT coefficients by using an absurd velocity
    do i = 1, 10
        params(i)%dopplerVelocityInput_old = -1e10
    enddo
    
    RETURN
end subroutine c_init

subroutine c_exit(index) bind(c)
integer(c_int), intent(in) :: index

    deallocate(observation(index)%wl)
    deallocate(inversion(index)%stokes_unperturbed)
    deallocate(fixed(index)%stokes_boundary)
    deallocate(fixed(index)%epsI)
    deallocate(fixed(index)%epsQ)
    deallocate(fixed(index)%epsU)
    deallocate(fixed(index)%epsV)
    deallocate(fixed(index)%etaI)
    deallocate(fixed(index)%etaQ)
    deallocate(fixed(index)%etaU)
    deallocate(fixed(index)%etaV)
    deallocate(fixed(index)%rhoQ)
    deallocate(fixed(index)%rhoU)
    deallocate(fixed(index)%rhoV)           
    deallocate(fixed(index)%dtau)

    deallocate(fixed(index)%epsilon)
    deallocate(fixed(index)%epsilon_zeeman)
    deallocate(fixed(index)%eta)
    deallocate(fixed(index)%eta_zeeman)
    deallocate(fixed(index)%eta_stim)
    deallocate(fixed(index)%eta_stim_zeeman)
    deallocate(fixed(index)%mag_opt)
    deallocate(fixed(index)%mag_opt_zeeman)
    deallocate(fixed(index)%mag_opt_stim)
    deallocate(fixed(index)%mag_opt_stim_zeeman)

    
end subroutine c_exit

end module pyHazelMod