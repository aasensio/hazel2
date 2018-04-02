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
    real(kind=8) :: StokesM(4), kappa_star(4,4), identity(4,4), source(4), m1(4,4), m2(4,4), Stokes0(4)
    real(kind=8) :: O_evol(4,4), psi_matrix(4,4), Stokes1(4)

    !
    ! JDLCR: init PSF?
    ! It will only be done the first time internally in the function
    !
        call init_psf()

        error = 0
        
! ! Fill the statistical equilibrium equations
!         call fill_SEE(in_params, in_fixed, 1, error)        

! ! If the solution of the SEE gives an error, return
!         if (error == 1) return
                
! ! Calculate the absorption/emission coefficients for a given transition
!         call calc_rt_coef(in_params, in_fixed, in_observation, 1)
                        
                
!****************       
! Slab case with EXACT SOLUTION
!****************
        if (.not.associated(in_fixed%epsI)) allocate(in_fixed%epsI(in_fixed%no))
        if (.not.associated(in_fixed%epsQ)) allocate(in_fixed%epsQ(in_fixed%no))
        if (.not.associated(in_fixed%epsU)) allocate(in_fixed%epsU(in_fixed%no))
        if (.not.associated(in_fixed%epsV)) allocate(in_fixed%epsV(in_fixed%no))
        if (.not.associated(in_fixed%etaI)) allocate(in_fixed%etaI(in_fixed%no))
        if (.not.associated(in_fixed%etaQ)) allocate(in_fixed%etaQ(in_fixed%no))
        if (.not.associated(in_fixed%etaU)) allocate(in_fixed%etaU(in_fixed%no))
        if (.not.associated(in_fixed%etaV)) allocate(in_fixed%etaV(in_fixed%no))
        if (.not.associated(in_fixed%rhoQ)) allocate(in_fixed%rhoQ(in_fixed%no))
        if (.not.associated(in_fixed%rhoU)) allocate(in_fixed%rhoU(in_fixed%no))
        if (.not.associated(in_fixed%rhoV)) allocate(in_fixed%rhoV(in_fixed%no))           
        if (.not.associated(in_fixed%dtau)) allocate(in_fixed%dtau(in_fixed%no))
        
        
        identity = 0.d0
        do i = 1, 4
            identity(i,i) = 1.d0
        enddo
        
        ! StokesM(1) = in_fixed%Stokes_incident(0)
        ! StokesM(2) = in_fixed%Stokes_incident(1)
        ! StokesM(3) = in_fixed%Stokes_incident(2)
        ! StokesM(4) = in_fixed%Stokes_incident(3)

                    
! Emission              
        in_fixed%epsI = in_fixed%epsilon(0,:)
        in_fixed%epsQ = in_fixed%epsilon(1,:)
        in_fixed%epsU = in_fixed%epsilon(2,:)
        in_fixed%epsV = in_fixed%epsilon(3,:)
        
! Absorption including stimulated emission
        in_fixed%etaI = in_fixed%eta(0,:) - use_stim_emission_RT * in_fixed%eta_stim(0,:)
        in_fixed%etaQ = in_fixed%eta(1,:) - use_stim_emission_RT * in_fixed%eta_stim(1,:)
        in_fixed%etaU = in_fixed%eta(2,:) - use_stim_emission_RT * in_fixed%eta_stim(2,:)
        in_fixed%etaV = in_fixed%eta(3,:) - use_stim_emission_RT * in_fixed%eta_stim(3,:)

! Magneto-optical effects
        if (use_mag_opt_RT == 1) then
            in_fixed%rhoQ = in_fixed%mag_opt(1,:) - use_stim_emission_RT * in_fixed%mag_opt_stim(1,:)
            in_fixed%rhoU = in_fixed%mag_opt(2,:) - use_stim_emission_RT * in_fixed%mag_opt_stim(2,:)
            in_fixed%rhoV = in_fixed%mag_opt(3,:) - use_stim_emission_RT * in_fixed%mag_opt_stim(3,:)
        else
            in_fixed%rhoQ = 0.d0
            in_fixed%rhoU = 0.d0
            in_fixed%rhoV = 0.d0
        endif

! Second component
        ds = in_params%dtau / maxval(in_fixed%etaI)
        in_fixed%dtau = in_fixed%etaI * ds
                
        do i = 1, in_fixed%no

            StokesM(1:4) = in_fixed%stokes_boundary(0:3,i)
            
            call fill_absorption_matrix(kappa_star,in_fixed%etaI(i),in_fixed%etaQ(i),in_fixed%etaU(i),in_fixed%etaV(i),in_fixed%rhoQ(i),in_fixed%rhoU(i),in_fixed%rhoV(i))
            kappa_star = kappa_star / in_fixed%etaI(i)
            source(1) = in_fixed%epsI(i) / in_fixed%etaI(i)
            source(2) = in_fixed%epsQ(i) / in_fixed%etaI(i)
            source(3) = in_fixed%epsU(i) / in_fixed%etaI(i)
            source(4) = in_fixed%epsV(i) / in_fixed%etaI(i)

! Evaluate the evolution operator
            call evol_operator(kappa_star,in_fixed%dtau(i),O_evol)

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

                                                
        ! if (in_observation%normalization == 'peak') then
        !     Imax = maxval(output(0,:))
        ! else
        !     Imax = output(0,1)
        ! endif
        ! do i = 0, 3
        !     output(i,:) = output(i,:) !/ Imax                
        ! enddo
            
        in_fixed%total_forward_modeling = in_fixed%total_forward_modeling + 1

        !
        ! JDLCR: convolve all Stokes Parameters with the spectral PSF.
        !
        call convolve(in_fixed%no, output)

    
    end subroutine do_synthesis
end module synth
