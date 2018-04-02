module synth
use vars
use SEE
use rt_coef
implicit none
  !
  ! JDLCR: new vars for the convolution with the PSF
  !
  integer :: npsf
  real(kind=8), allocatable :: psf(:)
  !
  PRIVATE :: npsf, psf
contains

  ! -------------------------------------------------------------------------
  !
  ! JDLCR: This function will only read the PSF in the first call.
  !        Right now hardwired to psf.txt.
  !
  ! -------------------------------------------------------------------------
  subroutine init_psf()
    implicit none
    character(len=7), parameter :: filename = 'psf.txt'
    integer :: unit, ii
    logical :: psf_exists

    if(allocated(psf)) return

    psf_exists = .FALSE.
    INQUIRE( FILE=filename, EXIST=psf_exists) 
    if(.not. psf_exists) return


    !
    ! Open PSF file and read
    !
    print *, 'Using PSF'
    unit = 1
    OPEN(unit, FILE=filename, status='OLD')
    read(unit,*) npsf ! Number of elements of the PSF


    allocate(psf(npsf)) ! allocate array to store the PSF

    
    do ii=1,npsf
       read(unit,*) psf(ii)
    end do

    
    CLOSE(unit)

  end subroutine init_psf

  ! -------------------------------------------------------------------------
  !
  ! JDLCR: This function will compute the convolution without FFTs.
  !        It will not pad the arrays, just use the part of the PSF that is
  !        inside the range of the spectra.
  !
  ! -------------------------------------------------------------------------
  subroutine convolve(n, sp)
    implicit none
    integer :: n, ii, jj, w0, w1, ww, npsf2, ss
    real(8) :: sp(0:3, n), res(n), psfsum


    if(.not. allocated(psf)) return

    npsf2 = npsf/2

    do ss = 0,3
       do ww = 1, n
          w0 = max(ww - npsf2, 1)
          w1 = min(ww + npsf2, n)
          ii = w0 - (ww - npsf2)
          jj = (ww + npsf2) - w1       
          res(ww) = sum(psf(1+ii:npsf-jj)*sp(ss,w0:w1)) / sum(psf(1+ii:npsf-jj))
       end do
       sp(ss,:) = res(:)
    end do

  end subroutine convolve

!------------------------------------------------------------
! Do a synthesis calling the appropriate routines
!------------------------------------------------------------
    subroutine do_synthesis(in_params,in_fixed,in_observation,output, error)
    type(variable_parameters) :: in_params, in_trial
    type(type_observation) :: in_observation
    type(fixed_parameters) :: in_fixed
    integer :: i, error
    real(kind=8) :: output(0:3,in_fixed%no), I0, Q0, U0, V0, ds, Imax, mu, Ic, factor, eta0, psim, psi0, sh, ds2
    real(kind=8) :: wstep
    real(kind=8), allocatable, dimension(:) :: epsI, epsQ, epsU, epsV, etaI, etaQ, etaU, etaV, dtau
    real(kind=8), allocatable, dimension(:) :: epsI2, epsQ2, epsU2, epsV2, etaI2, etaQ2, etaU2, etaV2
    real(kind=8), allocatable, dimension(:) :: rhoQ, rhoU, rhoV, delta
    real(kind=8), allocatable, dimension(:) :: rhoQ2, rhoU2, rhoV2
    real(kind=8), allocatable :: StokesM(:), kappa_prime(:,:), kappa_star(:,:), identity(:,:), source(:), m1(:,:), m2(:,:), Stokes0(:)
    real(kind=8), allocatable :: O_evol(:,:), psi_matrix(:,:), Stokes1(:)

    !
    ! JDLCR: init PSF?
    ! It will only be done the first time internally in the function
    !
        call init_psf()

        error = 0
        
! Fill the statistical equilibrium equations
        call fill_SEE(in_params, in_fixed, 1, error)        

! If the solution of the SEE gives an error, return
        if (error == 1) return
                
! Calculate the absorption/emission coefficients for a given transition
        call calc_rt_coef(in_params, in_fixed, in_observation, 1)
                        
                
!****************       
! Slab case with EXACT SOLUTION
!****************
        if (synthesis_mode == 5) then
            if (.not.allocated(epsI)) allocate(epsI(in_fixed%no))
            if (.not.allocated(epsQ)) allocate(epsQ(in_fixed%no))
            if (.not.allocated(epsU)) allocate(epsU(in_fixed%no))
            if (.not.allocated(epsV)) allocate(epsV(in_fixed%no))
            if (.not.allocated(etaI)) allocate(etaI(in_fixed%no))
            if (.not.allocated(etaQ)) allocate(etaQ(in_fixed%no))
            if (.not.allocated(etaU)) allocate(etaU(in_fixed%no))
            if (.not.allocated(etaV)) allocate(etaV(in_fixed%no))
            if (.not.allocated(rhoQ)) allocate(rhoQ(in_fixed%no))
            if (.not.allocated(rhoU)) allocate(rhoU(in_fixed%no))
            if (.not.allocated(rhoV)) allocate(rhoV(in_fixed%no))           
            if (.not.allocated(dtau)) allocate(dtau(in_fixed%no))
            
            if (.not.allocated(StokesM)) allocate(StokesM(4))
            if (.not.allocated(identity)) then
                allocate(identity(4,4))
                identity = 0.d0
                do i = 1, 4
                    identity(i,i) = 1.d0
                enddo
            endif
            if (.not.allocated(source)) allocate(source(4))
            if (.not.allocated(kappa_star)) allocate(kappa_star(4,4))
            if (.not.allocated(O_evol)) allocate(O_evol(4,4))
            if (.not.allocated(psi_matrix)) allocate(psi_matrix(4,4))
            if (.not.allocated(m1)) allocate(m1(4,4))
            if (.not.allocated(m2)) allocate(m2(4,4))
            if (.not.allocated(Stokes0)) allocate(Stokes0(4))
            if (.not.allocated(Stokes1)) allocate(Stokes1(4))
            
            ! StokesM(1) = in_fixed%Stokes_incident(0)
            ! StokesM(2) = in_fixed%Stokes_incident(1)
            ! StokesM(3) = in_fixed%Stokes_incident(2)
            ! StokesM(4) = in_fixed%Stokes_incident(3)

                        
            if (in_fixed%use_atomic_pol == 1 .or. in_fixed%use_atomic_pol == -1) then
! Emission              
                epsI = epsilon(0,:)
                epsQ = epsilon(1,:)
                epsU = epsilon(2,:)
                epsV = epsilon(3,:)
                
! Absorption including stimulated emission
                etaI = eta(0,:) - use_stim_emission_RT * eta_stim(0,:)
                etaQ = eta(1,:) - use_stim_emission_RT * eta_stim(1,:)
                etaU = eta(2,:) - use_stim_emission_RT * eta_stim(2,:)
                etaV = eta(3,:) - use_stim_emission_RT * eta_stim(3,:)

! Magneto-optical effects
                if (use_mag_opt_RT == 1) then
                    rhoQ = mag_opt(1,:) - use_stim_emission_RT * mag_opt_stim(1,:)
                    rhoU = mag_opt(2,:) - use_stim_emission_RT * mag_opt_stim(2,:)
                    rhoV = mag_opt(3,:) - use_stim_emission_RT * mag_opt_stim(3,:)
                else
                    rhoQ = 0.d0
                    rhoU = 0.d0
                    rhoV = 0.d0
                endif
            else
! Emission
                epsI = epsilon_zeeman(0,:)
                epsQ = epsilon_zeeman(1,:)
                epsU = epsilon_zeeman(2,:)
                epsV = epsilon_zeeman(3,:)

! Absorption including stimulated emission
                etaI = eta_zeeman(0,:) - use_stim_emission_RT * eta_stim_zeeman(0,:) + 1.d-20
                etaQ = eta_zeeman(1,:) - use_stim_emission_RT * eta_stim_zeeman(1,:) + 1.d-20
                etaU = eta_zeeman(2,:) - use_stim_emission_RT * eta_stim_zeeman(2,:) + 1.d-20
                etaV = eta_zeeman(3,:) - use_stim_emission_RT * eta_stim_zeeman(3,:) + 1.d-20

! Magneto-optical terms
                if (use_mag_opt_RT == 1) then
                    rhoQ = mag_opt_zeeman(1,:) - use_stim_emission_RT * mag_opt_stim_zeeman(1,:)
                    rhoU = mag_opt_zeeman(2,:) - use_stim_emission_RT * mag_opt_stim_zeeman(2,:)
                    rhoV = mag_opt_zeeman(3,:) - use_stim_emission_RT * mag_opt_stim_zeeman(3,:)
                else
                    rhoQ = 0.d0
                    rhoU = 0.d0
                    rhoV = 0.d0
                endif
            endif

! Second component
            ds = in_params%dtau / maxval(etaI)
            dtau = etaI * ds
            
            if (in_params%nslabs == 2) then
! Now the Doppler shift
! First calculate the difference in wavelength between the two components
                sh = in_fixed%wl * (in_params%vmacro2 - in_params%vmacro) * 1.d5 / PC

! Then calculate the number of "pixels" to carry out the correct Doppler shift

! If inverting, assume that the step in wavelength is constant and obtain it from
! the observations
                if (working_mode == 1) then
                    wstep = -(in_observation%wl(2)-in_observation%wl(1))
                endif

! If synthesizing, calculate it because the step is constant
                if (working_mode == 0) then
                    wstep = (in_fixed%omax - in_fixed%omin) / in_fixed%no
                endif
                
                sh = sh / wstep

                ds2 = in_params%dtau2 / maxval(etaI)
                dtau = etaI * ds + fft_shift(etaI, sh) * ds2

                epsI = in_params%dtau * epsI + in_params%dtau2 * fft_shift(epsI, sh)
                epsQ = in_params%dtau * epsQ + in_params%dtau2 * fft_shift(epsQ, sh)
                epsU = in_params%dtau * epsU + in_params%dtau2 * fft_shift(epsU, sh)
                epsV = in_params%dtau * epsV + in_params%dtau2 * fft_shift(epsV, sh)
                etaI = in_params%dtau * etaI + in_params%dtau2 * fft_shift(etaI, sh)
                etaQ = in_params%dtau * etaQ + in_params%dtau2 * fft_shift(etaQ, sh)
                etaU = in_params%dtau * etaU + in_params%dtau2 * fft_shift(etaU, sh)
                etaV = in_params%dtau * etaV + in_params%dtau2 * fft_shift(etaV, sh)
                rhoQ = in_params%dtau * rhoQ + in_params%dtau2 * fft_shift(rhoQ, sh)
                rhoU = in_params%dtau * rhoU + in_params%dtau2 * fft_shift(rhoU, sh)
                rhoV = in_params%dtau * rhoV + in_params%dtau2 * fft_shift(rhoV, sh)

            endif           
            
            do i = 1, in_fixed%no

                StokesM(1:4) = in_fixed%stokes_boundary(0:3,i)
                
                call fill_absorption_matrix(kappa_star,etaI(i),etaQ(i),etaU(i),etaV(i),rhoQ(i),rhoU(i),rhoV(i))
                kappa_star = kappa_star / etaI(i)
                source(1) = epsI(i) / etaI(i)
                source(2) = epsQ(i) / etaI(i)
                source(3) = epsU(i) / etaI(i)
                source(4) = epsV(i) / etaI(i)

! Evaluate the evolution operator
                call evol_operator(kappa_star,dtau(i),O_evol)

! Calculate K*^(-1)
                m1 = kappa_star
                call invert(m1)

                m2 = identity - O_evol
                Psi_matrix = matmul(m1,m2)

! Simplified version taking into account that the source function is constant, so that
!  I0 = exp(-K^* * tau_MO) * I_sun + (PsiM+Psi0)*S  with
! PsiM = U0-U1/tau_MO    and PsiO = U1/m
! U0 = (K*)^(-1) (1-exp(-K^* tau_MO)    and      U1 = (K*)^(-1) (m*1 - U0)
                Stokes0 = matmul(O_evol,StokesM) + matmul(Psi_matrix,source * in_params%beta)
                
                output(0,i) = Stokes0(1)
                output(1,i) = Stokes0(2)
                output(2,i) = Stokes0(3)
                output(3,i) = Stokes0(4)
                
            enddo

! If we include an additional slab with different field
            if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then

! Fill the statistical equilibrium equations
                call fill_SEE(in_params, in_fixed, 2, error)

! If the solution of the SEE gives an error, return
                if (error == 1) return
                
! Calculate the absorption/emission coefficients for a given transition
                call calc_rt_coef(in_params, in_fixed, in_observation, 2)

                if (in_fixed%use_atomic_pol == 1 .or. in_fixed%use_atomic_pol == -1) then
    ! Emission              
                    epsI = epsilon(0,:)
                    epsQ = epsilon(1,:)
                    epsU = epsilon(2,:)
                    epsV = epsilon(3,:)
                    
    ! Absorption including stimulated emission
                    etaI = eta(0,:) - use_stim_emission_RT * eta_stim(0,:)
                    etaQ = eta(1,:) - use_stim_emission_RT * eta_stim(1,:)
                    etaU = eta(2,:) - use_stim_emission_RT * eta_stim(2,:)
                    etaV = eta(3,:) - use_stim_emission_RT * eta_stim(3,:)

    ! Magneto-optical effects
                    if (use_mag_opt_RT == 1) then
                        rhoQ = mag_opt(1,:) - use_stim_emission_RT * mag_opt_stim(1,:)
                        rhoU = mag_opt(2,:) - use_stim_emission_RT * mag_opt_stim(2,:)
                        rhoV = mag_opt(3,:) - use_stim_emission_RT * mag_opt_stim(3,:)
                    else
                        rhoQ = 0.d0
                        rhoU = 0.d0
                        rhoV = 0.d0
                    endif
                else
! Emission
                    epsI = epsilon_zeeman(0,:)
                    epsQ = epsilon_zeeman(1,:)
                    epsU = epsilon_zeeman(2,:)
                    epsV = epsilon_zeeman(3,:)

! Absorption including stimulated emission
                    etaI = eta_zeeman(0,:) - use_stim_emission_RT * eta_stim_zeeman(0,:) + 1.d-20
                    etaQ = eta_zeeman(1,:) - use_stim_emission_RT * eta_stim_zeeman(1,:) + 1.d-20
                    etaU = eta_zeeman(2,:) - use_stim_emission_RT * eta_stim_zeeman(2,:) + 1.d-20
                    etaV = eta_zeeman(3,:) - use_stim_emission_RT * eta_stim_zeeman(3,:) + 1.d-20

! Magneto-optical terms
                    if (use_mag_opt_RT == 1) then
                        rhoQ = mag_opt_zeeman(1,:) - use_stim_emission_RT * mag_opt_stim_zeeman(1,:)
                        rhoU = mag_opt_zeeman(2,:) - use_stim_emission_RT * mag_opt_stim_zeeman(2,:)
                        rhoV = mag_opt_zeeman(3,:) - use_stim_emission_RT * mag_opt_stim_zeeman(3,:)
                    else
                        rhoQ = 0.d0
                        rhoU = 0.d0
                        rhoV = 0.d0
                    endif
                endif

! Second component
                ds = in_params%dtau2 / maxval(etaI)
                dtau = etaI * ds

                do i = 1, in_fixed%no
                    call fill_absorption_matrix(kappa_star,etaI(i),etaQ(i),etaU(i),etaV(i),rhoQ(i),rhoU(i),rhoV(i))
                    kappa_star = kappa_star / etaI(i)
                    source(1) = epsI(i) / etaI(i)
                    source(2) = epsQ(i) / etaI(i)
                    source(3) = epsU(i) / etaI(i)
                    source(4) = epsV(i) / etaI(i)

    ! Evaluate the evolution operator
                    call evol_operator(kappa_star,dtau(i),O_evol)

    ! Calculate K*^(-1)
                    m1 = kappa_star
                    call invert(m1)

                    m2 = identity - O_evol
                    Psi_matrix = matmul(m1,m2)

! Simplified version taking into account that the source function is constant, so that
!  I0 = exp(-K^* * tau_MO) * I_sun + (PsiM+Psi0)*S  with
! PsiM = U0-U1/tau_MO    and PsiO = U1/m
! U0 = (K*)^(-1) (1-exp(-K^* tau_MO)    and      U1 = (K*)^(-1) (m*1 - U0)

! Two components one after the other
                    if (in_params%nslabs == 3) then
                        StokesM(1:4) = output(0:3,i)

                        Stokes0 = matmul(O_evol,StokesM) + matmul(Psi_matrix,source * in_params%beta2)

                        output(0,i) = Stokes0(1)
                        output(1,i) = Stokes0(2)
                        output(2,i) = Stokes0(3)
                        output(3,i) = Stokes0(4)
                    endif

! Two components side by side inside the pixel with a filling factor
                    if (in_params%nslabs == -2) then

                        Stokes1 = matmul(O_evol,StokesM) + matmul(Psi_matrix,source * in_params%beta2)

                        output(0,i) = in_params%ff * output(0,i) + (1.d0-in_params%ff) * Stokes1(1)
                        output(1,i) = in_params%ff * output(1,i) + (1.d0-in_params%ff) * Stokes1(2)
                        output(2,i) = in_params%ff * output(2,i) + (1.d0-in_params%ff) * Stokes1(3)
                        output(3,i) = in_params%ff * output(3,i) + (1.d0-in_params%ff) * Stokes1(4)
                    endif
                
                enddo

            endif
                                                    
            if (in_observation%normalization == 'peak') then
                Imax = maxval(output(0,:))
            else
                Imax = output(0,1)
            endif
            do i = 0, 3
                output(i,:) = output(i,:) !/ Imax                
            enddo
            
        endif

        in_fixed%total_forward_modeling = in_fixed%total_forward_modeling + 1

        !
        ! JDLCR: convolve all Stokes Parameters with the spectral PSF.
        !
        call convolve(in_fixed%no, output)

    
    end subroutine do_synthesis
end module synth
