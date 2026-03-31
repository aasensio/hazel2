module nonlinear_transfer
use vars
use SEE
use rt_coef
use maths
implicit none
contains

!------------------------------------------------------------
! Read the file with the slab data
!------------------------------------------------------------
	subroutine set_slab(slab, logn, dz, v, B, thB, chiB)	
	type(type_slab) :: slab    
    real(kind=8) :: logn, dz, dz_slice, v, B, thB, chiB
	integer :: i
			
        dz_slice = dz / real(slab%nshells)
        do i = 1, slab%nshells
            slab%z(i) = dz_slice * real(i - 1.0) * 1.d5
            slab%density(i) = 10.0**logn
        enddo

        slab%velocity = v
        slab%B = B
        slab%thB = thB
        slab%chiB = chiB
        
	end subroutine set_slab

!------------------------------------------------------------
! Do a synthesis calling the appropriate routines
!------------------------------------------------------------
    subroutine do_transfer(in_params, in_fixed, in_observation, output, error)
    type(variable_parameters) :: in_params, in_trial
    type(type_observation) :: in_observation
    type(fixed_parameters) :: in_fixed
    type(type_slab) :: slab
    integer :: i, loop_iteration, loop_shell, error, j
    real(kind=8) :: output(0:3,in_fixed%no), I0, Q0, U0, V0, ds, Imax, mu, Ic, factor, eta0, psim, psi0
    real(kind=8) :: illumination_cone_cosine, illumination_cone_sine, dnum, relative_change(2), vmacro
    real(kind=8), allocatable, dimension(:) :: epsI, epsQ, epsU, epsV, etaI, etaQ, etaU, etaV, dtau            
    real(kind=8), allocatable, dimension(:) :: rhoQ, rhoU, rhoV, delta, prof(:)
    real(kind=8), allocatable :: StokesM(:), kappa_prime(:,:), kappa_star(:,:), identity(:,:)
    real(kind=8), allocatable :: source(:), m1(:,:), m2(:,:), Stokes0(:)
    real(kind=8), allocatable :: O_evol(:,:), psi_matrix(:,:), J00(:), J20(:), J00_nu(:,:), J20_nu(:,:)

        error = 0

        slab%nshells = in_fixed%nslabs
		
		slab%nmus_photosphere = in_fixed%nmus_photo
        slab%nmus_nophotosphere = in_fixed%nmus_nophoto
        		        
		slab%nmus = slab%nmus_photosphere + slab%nmus_nophotosphere

		allocate(slab%mus(slab%nmus))
		allocate(slab%weights(slab%nmus))
		allocate(slab%z(slab%nshells))
		allocate(slab%density(slab%nshells))
		allocate(slab%velocity(slab%nshells))
		allocate(slab%B(slab%nshells))
		allocate(slab%thB(slab%nshells))
		allocate(slab%chiB(slab%nshells))

        call set_slab(slab, in_params%logn, in_fixed%dz, in_params%vmacro, in_params%bgauss, in_params%thetabd, in_params%chibd)

        print *, 'Starting transfer...'
        
        allocate(J00_nu(slab%nshells,in_fixed%no))
        allocate(J20_nu(slab%nshells,in_fixed%no))
        allocate(J00(slab%nshells))
        allocate(J20(slab%nshells))
        
! Calculate the illumination cone at a given height for calculating the
! angular integration
        illumination_cone_sine = RSUN / (RSUN + in_params%height)
        illumination_cone_cosine = sqrt(1.d0-illumination_cone_sine**2)

! Set weights and mus for the photospheric cone
        call gauleg(illumination_cone_cosine,1.d0,slab%mus(1:slab%nmus_photosphere),&
            slab%weights(1:slab%nmus_photosphere),slab%nmus_photosphere)
        
! Set weights and mus for the rest of angles
        call gauleg(-1.d0,illumination_cone_cosine,slab%mus(slab%nmus_photosphere+1:slab%nmus),&
            slab%weights(slab%nmus_photosphere+1:slab%nmus),slab%nmus_nophotosphere)

        allocate(slab%nbar(slab%nshells,atom%ntran))
        allocate(slab%omega(slab%nshells,atom%ntran))

        allocate(slab%nbar_old(slab%nshells,atom%ntran))
        allocate(slab%omega_old(slab%nshells,atom%ntran))

        allocate(slab%emission_vector(4,slab%nshells,in_fixed%no))
        allocate(slab%propagation_matrix(4,4,slab%nshells,in_fixed%no))
        allocate(slab%boundary(4,slab%nmus,in_fixed%no))

        allocate(slab%tau(slab%nshells,in_fixed%no))

        allocate(prof(in_fixed%no))

        do i = 1, atom%ntran
            
! Initialize the values of nbar and omega to the Allen's values to start iterating
! Use Allen's tables to calculate the anisotropy and the value of nbar
            ! slab%nbar(:,i) = nbar_allen(atom%wavelength(i), in_fixed, in_params, atom%reduction_factor(i) * in_fixed%nbarExternal(i))
            ! slab%omega(:,i) = omega_allen(atom%wavelength(i), in_fixed, in_params, atom%reduction_factor_omega(i) * in_fixed%omegaExternal(i))

            slab%nbar(:,i) = nbar_allen(atom%wavelength(i), in_fixed, in_params, 1.d0)
            slab%omega(:,i) = omega_allen(atom%wavelength(i), in_fixed, in_params, 1.d0)

        enddo

! Boundary conditions
        slab%boundary = 0.d0
        do i = 1, slab%nmus_photosphere
! From Eq. 12.33 of Landi degl'Innocenti & Landolfi (2004)
            slab%boundary(1,i,:) = I0_allen(in_fixed%wl, sqrt(slab%mus(i)**2-illumination_cone_cosine**2) / illumination_cone_sine)
        enddo

        slab%nbar_old = slab%nbar
        slab%omega_old = slab%omega

        relative_change = 100.d0

        loop_iteration = 1

! velocity-free approximation
        vmacro = in_params%vmacro
        in_params%vmacro = 0.0
        
        do while (loop_iteration < 50 .and. relative_change(1) > 0.001d0 .and. relative_change(2) > 0.001d0)

            do loop_shell = 1, slab%nshells
            
! Fill and solve the statistical equilibrium equations
                in_params%bgauss = slab%B(i)
                in_params%thetabd = slab%thB(i)
                in_params%chibd = slab%chiB(i)

                nbarExternal = slab%nbar(loop_shell,:)
                omegaExternal = slab%omega(loop_shell,:)

                call fill_SEE(in_params, in_fixed, 1, .TRUE.)
                
! Calculate the absorption/emission coefficients for a given transition
! TODO : set velocity to zero for velocity-free approximation
                call calc_rt_coef(in_params, in_fixed, in_observation, 1, .TRUE.)

                                       
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

                if (.not.allocated(source)) allocate(source(4))
                if (.not.allocated(kappa_star)) allocate(kappa_star(4,4))

                StokesM(1) = in_fixed%Stokes_incident(0)
                StokesM(2) = in_fixed%Stokes_incident(1)
                StokesM(3) = in_fixed%Stokes_incident(2)
                StokesM(4) = in_fixed%Stokes_incident(3)
                
! Emission                
                epsI = in_fixed%epsilon(0,:)
                epsQ = in_fixed%epsilon(1,:)
                epsU = in_fixed%epsilon(2,:)
                epsV = in_fixed%epsilon(3,:)

            
! Absorption including stimulated emission
                etaI = in_fixed%eta(0,:) - in_fixed%eta_stim(0,:)
                etaQ = in_fixed%eta(1,:) - in_fixed%eta_stim(1,:)
                etaU = in_fixed%eta(2,:) - in_fixed%eta_stim(2,:)
                etaV = in_fixed%eta(3,:) - in_fixed%eta_stim(3,:)

! Magneto-optical effects                
                rhoQ = in_fixed%mag_opt(1,:) - in_fixed%mag_opt_stim(1,:)
                rhoU = in_fixed%mag_opt(2,:) - in_fixed%mag_opt_stim(2,:)
                rhoV = in_fixed%mag_opt(3,:) - in_fixed%mag_opt_stim(3,:)
                                
                ! Set emission vector at this shell
                slab%emission_vector(1,loop_shell,:) = epsI
                slab%emission_vector(2,loop_shell,:) = epsQ
                slab%emission_vector(3,loop_shell,:) = epsU
                slab%emission_vector(4,loop_shell,:) = epsV
                
! Set propagation matrix at this shell
                do i = 1, in_fixed%no
                    call fill_absorption_matrix(slab%propagation_matrix(:,:,loop_shell,i),&
                        etaI(i), etaQ(i), etaU(i), etaV(i), rhoQ(i), rhoU(i), rhoV(i))
                enddo
                
! Multiply by the density                
                slab%propagation_matrix(:,:,loop_shell,:) = slab%propagation_matrix(:,:,loop_shell,:) * slab%density(loop_shell)
                slab%emission_vector(:,loop_shell,:) = slab%emission_vector(:,loop_shell,:) * slab%density(loop_shell)                
                
            enddo
            
            J00 = 0.d0
            J20 = 0.d0
            J00_nu = 0.d0
            J20_nu = 0.d0

! Generate line profile
            dnum = in_params%vdopp*1.d5 / (in_fixed%wl*1.d-8)
            prof = profile(in_params%damping,in_observation%freq / dnum) / (dnum*SQRTPI)

! Solve the RT equation and calculate the tensors J00 and J20
            do i = 1, in_fixed%no
                call calculate_tensors(slab, i, J00_nu(:,i), J20_nu(:,i), in_fixed%gammad)
            enddo
            
! Carry out the integration over frequency weighted by the line profile                
            do i = 1, slab%nshells
                J00(i) = J00(i) + int_tabulated(-in_observation%freq, J00_nu(i,:)*prof)
                J20(i) = J20(i) + int_tabulated(-in_observation%freq, J20_nu(i,:)*prof)
            enddo
                                    
! Put the new values of nbar and omega
            slab%nbar(:,1) = J00 * (in_fixed%wl*1.d-8)**3 / (2.d0*PC*PH)
            slab%omega(:,1) = J20 / J00 * sqrt(2.d0)

            loop_iteration = loop_iteration + 1

            relative_change(1) = maxval(abs(slab%nbar - slab%nbar_old) / abs(slab%nbar))
            relative_change(2) = maxval(abs(slab%omega - slab%omega_old) / abs(slab%omega))

            slab%nbar_old = slab%nbar
            slab%omega_old = slab%omega

            print *, 'Maximum relative change : ', relative_change

        enddo

! Recompute the propagation matrix and emission vector for the last time
! for the synthesis of the Stokes profiles

        do loop_shell = 1, slab%nshells

! Now, restore the velocity for the final formal solution        
            in_params%vmacro = slab%velocity(loop_shell)
            in_params%bgauss = slab%B(loop_shell)
            in_params%thetabd = slab%thB(loop_shell)
            in_params%chibd = slab%chiB(loop_shell)            
            
! Fill and solve the statistical equilibrium equations
            nbarExternal = slab%nbar(loop_shell,:)
            omegaExternal = slab%omega(loop_shell,:)

            call fill_SEE(in_params, in_fixed, 1, .TRUE.)
                
! Calculate the absorption/emission coefficients for a given transition
! TODO : set velocity to zero for velocity-free approximation
            call calc_rt_coef(in_params, in_fixed, in_observation, 1, .TRUE.)
                                                
            if (.not.allocated(StokesM)) allocate(StokesM(4))

            if (.not.allocated(source)) allocate(source(4))
            if (.not.allocated(kappa_star)) allocate(kappa_star(4,4))

            StokesM(1) = in_fixed%Stokes_incident(0)
            StokesM(2) = in_fixed%Stokes_incident(1)
            StokesM(3) = in_fixed%Stokes_incident(2)
            StokesM(4) = in_fixed%Stokes_incident(3)
            
! Emission                
            epsI = in_fixed%epsilon(0,:)
            epsQ = in_fixed%epsilon(1,:)
            epsU = in_fixed%epsilon(2,:)
            epsV = in_fixed%epsilon(3,:)
            
! Absorption including stimulated emission
            etaI = in_fixed%eta(0,:) - in_fixed%eta_stim(0,:)
            etaQ = in_fixed%eta(1,:) - in_fixed%eta_stim(1,:)
            etaU = in_fixed%eta(2,:) - in_fixed%eta_stim(2,:)
            etaV = in_fixed%eta(3,:) - in_fixed%eta_stim(3,:)

            rhoQ = in_fixed%mag_opt(1,:) - in_fixed%mag_opt_stim(1,:)
            rhoU = in_fixed%mag_opt(2,:) - in_fixed%mag_opt_stim(2,:)
            rhoV = in_fixed%mag_opt(3,:) - in_fixed%mag_opt_stim(3,:)
            

! Set emission vector at this shell
            slab%emission_vector(1,loop_shell,:) = epsI
            slab%emission_vector(2,loop_shell,:) = epsQ
            slab%emission_vector(3,loop_shell,:) = epsU
            slab%emission_vector(4,loop_shell,:) = epsV

! Set propagation matrix at this shell
            do i = 1, in_fixed%no
                call fill_absorption_matrix(slab%propagation_matrix(:,:,loop_shell,i),&
                    etaI(i), etaQ(i), etaU(i), etaV(i), rhoQ(i), rhoU(i), rhoV(i))
            enddo

! Multiply by the density
            slab%propagation_matrix(:,:,loop_shell,:) = slab%propagation_matrix(:,:,loop_shell,:) * slab%density(loop_shell)
            slab%emission_vector(:,loop_shell,:) = slab%emission_vector(:,loop_shell,:) * slab%density(loop_shell)

        enddo

        call synthesize_stokes(slab, in_fixed, output)

        open(unit=18,file='Jbar_tensors.dat',action='write',status='replace')
        write(18,*) slab%nshells
        do i = 1, slab%nshells
            write(18,*) maxval(slab%tau(i,:)), slab%nbar(i,1), slab%omega(i,1)
        enddo
        close(18)

        deallocate(J00_nu)
        deallocate(J20_nu)
        deallocate(J00)
        deallocate(J20)
        deallocate(prof)

        deallocate(slab%mus)
		deallocate(slab%weights)
		deallocate(slab%z)
		deallocate(slab%density)
		deallocate(slab%velocity)
		deallocate(slab%B)
		deallocate(slab%thB)
		deallocate(slab%chiB)
            
    end subroutine do_transfer

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------    
    subroutine calculate_tensors(slab, freq, J00, J20, gamma)
    type(type_slab) :: slab
    integer :: freq
    real(kind=8) :: formal_sol_polarized(4), Inten(4), J00(:), J20(:), gamma
    integer :: k, km, kp, loop_mu, loop_direction
    real(kind=8) :: chim, chi0, chip, dtp, dtm, exu
    real(kind=8) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp, mu, Qtilde
    
    integer :: i, j, n, nmus, kfrom, kto, kstep
    real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:)
    real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4)
    
        n = slab%nshells
        nmus = slab%nmus

        ! 4x4 identity matrix
        identity_4x4 = 0.d0
        do i = 1, 4
            identity_4x4(i,i) = 1.d0
        enddo
        
        allocate(ab_matrix(4,4,n))
        allocate(source_vector(4,n))
                                
! Transform K into K* and then into K'
        do i = 1, 4
            do j = 1, 4
                ab_matrix(i,j,:) = slab%propagation_matrix(i,j,:,freq) / slab%propagation_matrix(1,1,:,freq)
            enddo
            ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
            source_vector(i,:) = slab%emission_vector(i,:,freq) / slab%propagation_matrix(1,1,:,freq)
        enddo
                
        J00 = 0.d0
        J20 = 0.d0

        do loop_mu = 1, nmus

            mu = slab%mus(loop_mu)
            
! If going upwards
            if (mu > 0) then
                kfrom = 2
                kto = n
                kstep = 1
            else
! If going downwards
                kfrom = n-1
                kto = 1
                kstep = -1
            endif

! Boundary condition
            Inten = slab%boundary(:,loop_mu,freq)            

! Calculate J00
            J00(kfrom-kstep) = J00(kfrom-kstep) + J00_FACTOR * slab%weights(loop_mu) * Inten(1)

! Calculate J20 taking into account the contribution of I, Q and U
            Qtilde = cos(2.d0*gamma*PI/180.d0) * Inten(2) - sin(2.d0*gamma*PI/180.d0) * Inten(3)
            J20(kfrom-kstep) = J20(kfrom-kstep) + J20_FACTOR * slab%weights(loop_mu) * &
                ( (3.d0*mu**2-1.d0) * Inten(1) - 3.d0*(1.d0-mu**2) * Qtilde )
            
            do k = kfrom, kto

! Parabolic short-characteristics
                if (k /= kto) then
                    km = k - 1
                    kp = k + 1
                    chim = slab%propagation_matrix(1,1,km,freq)                    
                    chi0 = slab%propagation_matrix(1,1,k,freq)                    
                    chip = slab%propagation_matrix(1,1,kp,freq)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)
                    chip = 0.d0
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = 0.d0
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = 0.d0
                endif
                    
                dtm = 0.5d0 * (chim + chi0) * dm
                dtp = 0.5d0 * (chi0 + chip) * dp
                                   
                if (dtm >= 1.d-4) then
                    exu = dexp(-dtm)
                else
                    exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
                endif
                
                call lin_sc(dtm,psim_lin,psi0_lin)
                mat1 = exu * identity_4x4 - psim_lin*ab_matrix(:,:,km)
                mat2 = identity_4x4 + psi0_lin * ab_matrix(:,:,k)
                call invert(mat2)
                        
                if (k /= kto) then
                    call par_sc(dtm,dtp,psim,psi0,psip)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
                else
                    call lin_sc(dtm,psim,psi0)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
                endif

! Calculate J00
                J00(k) = J00(k) + J00_FACTOR * slab%weights(loop_mu) * Inten(1)

! Calculate J20 taking into account the contribution of I, Q and U
                Qtilde = cos(2.d0*gamma*PI/180.d0) * Inten(2) - sin(2.d0*gamma*PI/180.d0) * Inten(3)
                J20(k) = J20(k) + J20_FACTOR * slab%weights(loop_mu) * &
                    ( (3.d0*mu**2-1.d0) * Inten(1) - 3.d0*(1.d0-mu**2) * Qtilde )

                ! print *, k, Inten(1)
                        
            enddo
        enddo

        deallocate(ab_matrix)
        deallocate(source_vector)

    end subroutine calculate_tensors

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------    
    subroutine synthesize_stokes(slab, in_fixed, output)
    type(type_slab) :: slab
    type(fixed_parameters) :: in_fixed
    real(kind=8) :: output(0:3,in_fixed%no)
    real(kind=8) :: formal_sol_polarized(4), Inten(4)
    integer :: k, km, kp, freq
    real(kind=8) :: chim, chi0, chip, dtp, dtm, exu
    real(kind=8) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp, mu
    
    integer :: i, j, n, nmus, kfrom, kto, kstep
    real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:), total_tau(:)
    real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4), Ic
    
        n = slab%nshells
        nmus = slab%nmus

        ! 4x4 identity matrix
        identity_4x4 = 0.d0
        do i = 1, 4
            identity_4x4(i,i) = 1.d0
        enddo
        
        allocate(ab_matrix(4,4,n))
        allocate(source_vector(4,n))
        allocate(total_tau(in_fixed%no))

        total_tau = 0.d0
        slab%tau = 0.d0

        do freq = 1, in_fixed%no
                                
! Transform K into K* and then into K'
            do i = 1, 4
                do j = 1, 4
                    ab_matrix(i,j,:) = slab%propagation_matrix(i,j,:,freq) / slab%propagation_matrix(1,1,:,freq)
                enddo
                ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
                source_vector(i,:) = slab%emission_vector(i,:,freq) / slab%propagation_matrix(1,1,:,freq)
            enddo

            mu = 1.d0

! If going upwards
            kfrom = 2
            kto = n
            kstep = 1

! Boundary condition
            Inten = 0.d0
            Inten(1) = I0_allen(in_fixed%wl, mu)

            do k = kfrom, kto

! Parabolic short-characteristics
                if (k /= kto) then
                    km = k - 1
                    kp = k + 1
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)                    
                    chip = slab%propagation_matrix(1,1,kp,freq)                    
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = source_vector(:,kp)
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = dabs((slab%z(kp) - slab%z(k)) / mu)
                else
! Linear short-characteristics            
                    km = k - 1
                    chim = slab%propagation_matrix(1,1,km,freq)
                    chi0 = slab%propagation_matrix(1,1,k,freq)
                    chip = 0.d0
                    sm = source_vector(:,km)
                    s0 = source_vector(:,k)
                    sp = 0.d0
                    dm = dabs((slab%z(k) - slab%z(km)) / mu)
                    dp = 0.d0
                endif
                    
                dtm = 0.5d0 * (chim + chi0) * dm
                dtp = 0.5d0 * (chi0 + chip) * dp

                total_tau(freq) = total_tau(freq) + dtm
                slab%tau(k,freq) = slab%tau(km,freq) + dtm
                    
                if (dtm >= 1.d-4) then
                    exu = dexp(-dtm)
                else
                    exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
                endif
                
                call lin_sc(dtm,psim_lin,psi0_lin)
                mat1 = exu * identity_4x4 - psim_lin*ab_matrix(:,:,km)
                mat2 = identity_4x4 + psi0_lin * ab_matrix(:,:,k)
                call invert(mat2)
                        
                if (k /= kto) then
                    call par_sc(dtm,dtp,psim,psi0,psip)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
                else
                    call lin_sc(dtm,psim,psi0)
                    Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
                endif
                        
            enddo

            output(0:3,freq) = Inten
            
        enddo

        print *, 'Maximum tau=', maxval(total_tau)
        open(unit=18,file='final_tau.dat',action='write',status='replace')
        write(18,*) maxval(total_tau)
        do i = 1, in_fixed%no
            write(18,*) total_tau(i)
        enddo
        close(18)

! Return I/Ic, Q/Ic, U/Ic and V/Ic
        ! Ic = output(0,1)
        ! do i = 0, 3
            ! output(i,:) = output(i,:) / Ic
        ! enddo

        deallocate(ab_matrix)
        deallocate(source_vector)

    end subroutine synthesize_stokes

end module nonlinear_transfer
