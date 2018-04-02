module mpi_routines
use vars
implicit none

contains

!------------------------------------------------------------
! Broadcast Allen's data
!------------------------------------------------------------
    subroutine bcast_allen
    integer :: ierr
    include 'mpif.h'

        call MPI_Bcast(allen_ic,43*2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(allen_cl,22*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_allen

!------------------------------------------------------------
! Broadcast main configuration
!------------------------------------------------------------
    subroutine bcast_main_conf
    integer :: ierr
    include 'mpif.h'
        
        call MPI_Bcast(verbose_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(linear_solver,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(synthesis_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(working_mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_main_conf

!------------------------------------------------------------
! Broadcast experiment
!------------------------------------------------------------
    subroutine bcast_experiment_conf
    integer :: ierr
    include 'mpif.h'

        call MPI_Bcast(isti,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(imag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(idep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%delta_collision,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(use_paschen_back,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        call MPI_Bcast(params%bgauss,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%thetabd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%chibd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_Bcast(params%height,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%nslabs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Bcast(params%bgauss2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(params%thetabd2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(params%chibd2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif

        call MPI_Bcast(params%dtau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (params%nslabs == 2 .or. params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Bcast(params%dtau2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
        if (params%nslabs == -2) then
            call MPI_Bcast(params%ff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
        
        call MPI_Bcast(params%beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%Stokes_incident,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%nemiss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%use_atomic_pol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%thetad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%chid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%gammad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%omin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%omax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%no,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%wl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%Stokes_incident_mode,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        
        call MPI_Bcast(params%vdopp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Bcast(params%vdopp2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
        
        call MPI_Bcast(params%damping,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_Bcast(params%vmacro,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        if (params%nslabs /= 1) then
            call MPI_Bcast(params%vmacro2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(params%beta2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        endif
        
        call MPI_Bcast(use_mag_opt_RT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(use_stim_emission_RT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%total_forward_modeling,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_experiment_conf


!------------------------------------------------------------
! Broadcast atom model
!------------------------------------------------------------
    subroutine bcast_atom_model(myrank)
    integer :: ierr, myrank
    include 'mpif.h'

        call MPI_Bcast(is2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(n_terms,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)     
        call MPI_Bcast(jlimit2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(nrhos,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%ntran,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (myrank /= 0) then
            allocate(lsto2(n_terms))
            allocate(energy(n_terms,0:jlimit2))
            allocate(ntab(nrhos))
            allocate(j2tab(nrhos))
            allocate(jp2tab(nrhos))
            allocate(ktab(nrhos))
            allocate(qtab(nrhos))
            allocate(irtab(nrhos))
            allocate(rnutab(nrhos))
            allocate(atom%nterml(atom%ntran))
            allocate(atom%ntermu(atom%ntran))
            allocate(atom%ae(atom%ntran))
            allocate(atom%wavelength(atom%ntran))
            allocate(atom%reduction_factor(atom%ntran))
            allocate(atom%reduction_factor_omega(atom%ntran))
            allocate(atom%j10(atom%ntran))
        endif

        call MPI_Bcast(lsto2,n_terms,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(energy(:,0:jlimit2),n_terms*(jlimit2+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_Bcast(ntab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(j2tab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(jp2tab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ktab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(qtab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(irtab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(rnutab,nrhos,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                
        call MPI_Bcast(atom%nterml,atom%ntran,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%ntermu,atom%ntran,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%ae,atom%ntran,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%wavelength,atom%ntran,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%reduction_factor,atom%ntran,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%reduction_factor_omega,atom%ntran,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atom%j10,atom%ntran,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(nbarExternal,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(omegaExternal,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_atom_model

!------------------------------------------------------------
! Broadcast main configuration
!------------------------------------------------------------
    subroutine bcast_parameters_to_invert(myrank)
    integer :: ierr, myrank
    include 'mpif.h'

        call MPI_Bcast(params%n_total,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(inversion%iter_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(inversion%n_cycles,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        if (myrank /= 0) then
            allocate(params%inverted(params%n_total))
            allocate(inversion%cycles(params%n_total,inversion%n_cycles))
            allocate(inversion%stokes_weights(0:3,inversion%n_cycles))
            allocate(inversion%algorithm(inversion%n_cycles))
        endif

        call MPI_Bcast(inversion%cycles,params%n_total*inversion%n_cycles,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(inversion%stokes_weights(0:3,:),4*inversion%n_cycles,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(inversion%algorithm,inversion%n_cycles,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        call MPI_Bcast(params%inverted,params%n_total,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(params%n_inverted,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_parameters_to_invert

!------------------------------------------------------------
! Broadcast information about sizes of observations
!------------------------------------------------------------
    subroutine bcast_observation_info
    integer :: ierr
    include 'mpif.h'

        call MPI_Bcast(observation%n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%no,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
    end subroutine bcast_observation_info

!------------------------------------------------------------
! Broadcast DIRECT configuration
!------------------------------------------------------------
    subroutine bcast_direct_info(myrank)
    integer :: ierr, myrank
    include 'mpif.h'

        if (myrank /= 0) then
            allocate(fixed%upper_direct(18))
            allocate(fixed%lower_direct(18))
        endif
        
        call MPI_Bcast(fixed%upper_direct,18,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%lower_direct,18,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%DIRmaxf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(fixed%volper,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine bcast_direct_info
!------------------------------------------------------------
! Send a new observation to a slave
!------------------------------------------------------------
    subroutine send_observation(observation, packagesize, slave, index_obs)
    type(type_observation) :: observation
    integer :: slave, packagesize, pos, index_obs
    character(len=packagesize) :: buffer
    integer :: ierr
    include 'mpif.h'

        pos = 0     
        call MPI_Pack(index_obs, 1, MPI_INTEGER, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(observation%wl, observation%n, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(observation%stokes(0:3,:), 4*observation%n, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)     
        call MPI_Pack(observation%sigma(0:3,:), 4*observation%n, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(fixed%Stokes_incident,4, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)        
        call MPI_Pack(fixed%thetad,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(fixed%gammad,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(params%height,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%bgauss,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)        
        call MPI_Pack(params%thetabd,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)       
        call MPI_Pack(params%chibd,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%bgauss2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)           
            call MPI_Pack(params%thetabd2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)          
            call MPI_Pack(params%chibd2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif

        call MPI_Pack(params%dtau,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        if (params%nslabs == 2 .or. params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%dtau2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif
        
        if (params%nslabs == -2) then
            call MPI_Pack(params%ff,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif

        call MPI_Pack(params%vdopp,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%vdopp2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif

        call MPI_Pack(params%damping,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

        call MPI_Pack(params%beta,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

        call MPI_Pack(params%vmacro,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        if (params%nslabs /= 1) then
            call MPI_Pack(params%vmacro2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)                       
            call MPI_Pack(params%beta2,1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif
        
        call MPI_Send(buffer, packagesize, MPI_PACKED, slave, 10, MPI_COMM_WORLD, ierr)
        
    end subroutine send_observation

!------------------------------------------------------------
! Receive a new observation from the master
!------------------------------------------------------------
    subroutine receive_observation(observation, packagesize, stop_flag, index_obs)
    type(type_observation) :: observation
    integer :: packagesize, pos, stop_flag, index_obs
    character(len=packagesize) :: buffer
    integer :: ierr
    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
            
        call MPI_Recv(buffer, packagesize, MPI_PACKED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

! If tag is 10 -> unpack the observation
        if (status(MPI_TAG) == 10) then
            stop_flag = 0   
            pos = 0
            call MPI_Unpack(buffer, packagesize, pos, index_obs, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, observation%wl, observation%n, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)        
            call MPI_Unpack(buffer, packagesize, pos, observation%stokes(0:3,:), 4*observation%n, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)       
            call MPI_Unpack(buffer, packagesize, pos, observation%sigma(0:3,:), 4*observation%n, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, fixed%Stokes_incident, 4, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, fixed%thetad, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, fixed%gammad, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, params%height, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

            call MPI_Unpack(buffer, packagesize, pos, params%bgauss, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, params%thetabd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, params%chibd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            
            if (params%nslabs == 3 .or. params%nslabs == -2) then
                call MPI_Unpack(buffer, packagesize, pos, params%bgauss2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, packagesize, pos, params%thetabd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, packagesize, pos, params%chibd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            endif

            call MPI_Unpack(buffer, packagesize, pos, params%dtau, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

            if (params%nslabs == 2 .or. params%nslabs == 3 .or. params%nslabs == -2) then
                call MPI_Unpack(buffer, packagesize, pos, params%dtau2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)              
            endif

            if (params%nslabs == -2) then
                call MPI_Unpack(buffer, packagesize, pos, params%ff, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            endif

            call MPI_Unpack(buffer, packagesize, pos, params%vdopp, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
            if (params%nslabs == 3 .or. params%nslabs == -2) then
                call MPI_Unpack(buffer, packagesize, pos, params%vdopp2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            endif

            call MPI_Unpack(buffer, packagesize, pos, params%damping, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

            call MPI_Unpack(buffer, packagesize, pos, params%beta, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

            call MPI_Unpack(buffer, packagesize, pos, params%vmacro, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            
            if (params%nslabs /= 1) then
                call MPI_Unpack(buffer, packagesize, pos, params%vmacro2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                call MPI_Unpack(buffer, packagesize, pos, params%beta2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            endif

            return
        endif

! This is the kill flag
        if (status(MPI_TAG) == 999) then
            stop_flag = 1
            return
        endif
        
    end subroutine receive_observation

!------------------------------------------------------------
! Send a converged model to the master
!------------------------------------------------------------
    subroutine send_model(params, error, in_observation, in_inversion, packagesize, myrank, index_obs)
    type(variable_parameters) :: params, error
    type(type_observation) :: in_observation
    type(type_inversion) :: in_inversion
    integer :: packagesize, pos, myrank, index_obs
    character :: buffer(packagesize)
    integer :: ierr
    include 'mpif.h'

! Pack the model
        pos = 0     
        call MPI_Pack(myrank, 1, MPI_INTEGER, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(index_obs, 1, MPI_INTEGER, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
! Model parameters
        call MPI_Pack(params%bgauss, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%bgauss, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%thetabd, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%thetabd, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%chibd, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%chibd, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%height, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%height, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%nslabs, 1, MPI_INTEGER, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%nslabs, 1, MPI_INTEGER, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%dtau, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%dtau, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

        if (params%nslabs == 2 .or. params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%dtau2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%dtau2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif

        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%bgauss2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%bgauss2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            
            call MPI_Pack(params%thetabd2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%thetabd2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            
            call MPI_Pack(params%chibd2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%chibd2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif
        
        call MPI_Pack(params%vdopp, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%vdopp, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Pack(params%vdopp2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%vdopp2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif
        
        call MPI_Pack(params%damping, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%damping, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        call MPI_Pack(params%vmacro, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%vmacro, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

        call MPI_Pack(params%beta, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(error%beta, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        
        if (params%nslabs /= 1) then
            call MPI_Pack(params%vmacro2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%vmacro2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)

            call MPI_Pack(params%beta2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%beta2, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif

        if (params%nslabs == -2) then
            call MPI_Pack(params%ff, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
            call MPI_Pack(error%ff, 1, MPI_DOUBLE_PRECISION, buffer, packagesize, pos, MPI_COMM_WORLD, ierr)
        endif
        
! And the synthetic profiles
        call MPI_Pack(in_inversion%stokes_unperturbed(0:3,:), 4*in_observation%n, MPI_DOUBLE_PRECISION, buffer, &
            packagesize, pos, MPI_COMM_WORLD, ierr)
                
! Send it to the master
        call MPI_Send(buffer, packagesize, MPI_PACKED, 0, 13, MPI_COMM_WORLD,ierr)

    end subroutine send_model

!------------------------------------------------------------
! Receive a converged model from a slave
!------------------------------------------------------------
    subroutine receive_model(params, error, in_observation, in_inversion, packagesize, which_slave, index_obs)
    type(variable_parameters) :: params, error
    type(type_observation) :: in_observation
    type(type_inversion) :: in_inversion
    integer :: packagesize, pos, index_obs
    character :: buffer(packagesize)
    integer :: ierr, which_slave
    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
    
! Receive data from slave
        call MPI_Recv(buffer, packagesize, MPI_PACKED, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
        
! Unpack the model
        pos = 0
        call MPI_Unpack(buffer, packagesize, pos, which_slave, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, index_obs, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%bgauss, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%bgauss, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%thetabd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%thetabd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%chibd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%chibd, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%height, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%height, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%nslabs, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%nslabs, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%dtau, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%dtau, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

        if (params%nslabs == 2 .or. params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Unpack(buffer, packagesize, pos, params%dtau2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%dtau2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        endif

        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Unpack(buffer, packagesize, pos, params%bgauss2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%bgauss2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            
            call MPI_Unpack(buffer, packagesize, pos, params%thetabd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%thetabd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            
            call MPI_Unpack(buffer, packagesize, pos, params%chibd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%chibd2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        endif
        
        call MPI_Unpack(buffer, packagesize, pos, params%vdopp, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%vdopp, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        if (params%nslabs == 3 .or. params%nslabs == -2) then
            call MPI_Unpack(buffer, packagesize, pos, params%vdopp2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%vdopp2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        endif
        
        call MPI_Unpack(buffer, packagesize, pos, params%damping, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%damping, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(buffer, packagesize, pos, params%vmacro, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%vmacro, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

        call MPI_Unpack(buffer, packagesize, pos, params%beta, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        call MPI_Unpack(buffer, packagesize, pos, error%beta, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

        if (params%nslabs /= 1) then
            call MPI_Unpack(buffer, packagesize, pos, params%vmacro2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%vmacro2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

            call MPI_Unpack(buffer, packagesize, pos, params%beta2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%beta2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
                        
        endif
        
        if (params%nslabs == -2) then
            call MPI_Unpack(buffer, packagesize, pos, params%ff, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_Unpack(buffer, packagesize, pos, error%ff, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        endif

! And the synthetic observations
        call MPI_Unpack(buffer, packagesize, pos, in_inversion%stokes_unperturbed(0:3,:), 4*in_observation%n, &
            MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        
    end subroutine receive_model
    

!------------------------------------------------------------
! Send the kill signal to all slaves
!------------------------------------------------------------
    subroutine kill_slave(slave)
    integer :: slave
    integer :: ierr, i
    include 'mpif.h'

! Send an empty message with the tag 999
        call MPI_Send(0, 0, MPI_INTEGER, slave, 999, MPI_COMM_WORLD, ierr)
        
    end subroutine kill_slave

end module mpi_routines
