program multiterm
use vars
use maths
use io
use SEE
use rt_coef
use marquardt
use synth
use allen
#if defined(mpi)
use mpi_routines
#else
use mpi_routines_fake
#endif
use inversion_mod
implicit none

    integer :: successful, n_procs_done, slave, status_obs
    character(len=120) :: temporal_file
    integer :: mpi_status, myrank, nprocs, ierr, package_size_model, package_size_obs, kill, index_obs, i, error
    integer, allocatable :: slave_active(:)
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    integer :: values(8)

#if defined(mpi)
    include 'mpif.h'
    call MPI_INIT(mpi_status)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpi_status)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpi_status)
    write(*,FMT='(A,I3,A,I3)') 'Node ', myrank, '/', nprocs
#else
    integer :: MPI_COMM_WORLD
    myrank = 0
    nprocs = 0
#endif

    nbarExternal = 1.d0
    omegaExternal = 1.d0
    
            
    allocate(slave_active(nprocs-1))
    if (myrank == 0) then
        slave_active = 1
!       open(unit=24,file='logfile',action='write',status='replace')
    else
        slave_active = 0
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
            
    temporal_file = 'temporal.prof'
    
! Initialize the random number generator
    call random_seed
    
! Initialize Allen's data
! Allocate memory for arrays
    allocate(allen_ic(43,2))
    allocate(allen_cl(22,3))

! Read Allen's data if master
    if (myrank == 0) then
        call read_allen_data
    endif

! Broadcast Allen's data
    call bcast_allen
    
    
! Fill the factorial array
    call factrl
                
! Read the configuration file   if master
    if (myrank == 0) then
        call read_config
    endif

! Broadcast relevant information from the main configuration file
! to all nodes
    call bcast_main_conf
  
        
! Read the file with the experiment
    if (myrank == 0) then
        call read_experiment(input_experiment, params, fixed)
    endif

! Broadcast relevant information from the experiment configuration file
! to all nodes
    call bcast_experiment_conf
    
! Read the atomic model 
    if (myrank == 0) then
        call read_model_file(input_model_file)
    endif

! Broadcast relevant information from the atom model configuration file
! to all nodes
    call bcast_atom_model(myrank)
        
    inversion%min_chisq = 1.d10
    successful = 0
    
!*********************************
!** SYNTHESIS MODE
!*********************************  
    if (working_mode == 0 .and. nprocs == 0) then
        if (verbose_mode == 1) then
            print *, 'Working in SYNTHESIS mode'
        endif
        observation%n = fixed%no
        allocate(observation%wl(observation%n))
        if (.not.associated(inversion%stokes_unperturbed)) allocate(inversion%stokes_unperturbed(0:3,fixed%no))

        call set_boundary_condition(fixed, inversion)

        call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
        call write_final_profiles(output_inverted_profiles,observation,inversion)
    endif
    
!*********************************
!** INVERSION MODE
!*********************************
    if (working_mode > 0 .and. nprocs >= 0) then
        if (verbose_mode == 1 .and. myrank == 0) then
            print *, 'Working in INVERSION mode'
        endif
        
! Read the parameters to be inverted
        if (myrank == 0) then
            call read_parameters_to_invert(input_inverted_parameters, params, inversion, fixed)
        endif

! Broadcast relevant information from which parameters to invert to all nodes
        call bcast_parameters_to_invert(myrank)

        if (myrank == 0) then
            call preread_observations(input_observed_profiles, fixed, observation)
        endif
        
! Broadcast relevant information about observation size
        call bcast_observation_info

        allocate(observation%wl(observation%n))
        allocate(observation%stokes(0:3,observation%n))
        allocate(observation%sigma(0:3,observation%n))

! Compute size of package for models and observations
        package_size_model = 4*sizeof(params%nslabs) + 2*9*sizeof(params%bgauss) + 4*sizeof(observation%wl)
            
        if (params%nslabs == 2) then
            package_size_model = package_size_model + 2*3*sizeof(params%bgauss)
        endif
        if (params%nslabs == 3) then
            package_size_model = package_size_model + 2*7*sizeof(params%bgauss)
        endif
        if (params%nslabs == -2) then
            package_size_model = package_size_model + 2*8*sizeof(params%bgauss)
        endif               

        package_size_obs = sizeof(n_procs_done) + sizeof(observation%wl) + sizeof(observation%stokes(0:3,:)) + &
            sizeof(observation%sigma(0:3,:)) + sizeof(fixed%Stokes_incident) + sizeof(params%height) + &
            sizeof(fixed%thetad) + sizeof(fixed%gammad) + 8*sizeof(params%bgauss)
             
        if (params%nslabs == 2) then
            package_size_obs = package_size_obs + 2*sizeof(params%bgauss)
        endif
        if (params%nslabs == 3) then
            package_size_obs = package_size_obs + 7*sizeof(params%bgauss)
        endif
        if (params%nslabs == -2) then
            package_size_obs = package_size_obs + 8*sizeof(params%bgauss)
        endif
                            
        if (.not.associated(inversion%stokes_unperturbed)) allocate(inversion%stokes_unperturbed(0:3,observation%n))
        if (.not.associated(inversion%stokes_perturbed)) allocate(inversion%stokes_perturbed(0:3,observation%n))
        if (.not.associated(inversion%dydx)) allocate(inversion%dydx(0:3,observation%n,params%n_total))

! Read DIRECT configuration file
        if (myrank == 0) then
            call read_direct_conf(direct_ranges, fixed)
        endif

! Broadcast DIRECT configuration
        call bcast_direct_info(myrank)

        n_procs_done = starting_pixel
        status_obs = 1

! Start master/slave process if more than two processes are found
        if (nprocs >= 2) then
        
! Initial distribution of works
            if (myrank == 0) then

! Loop over all process in the initial state distributing processes
                do while (status_obs == 1 .and. (n_procs_done - starting_pixel + 1) <= nprocs-1)                    
                    call read_observation(input_observed_profiles, fixed, observation, params, n_procs_done, status_obs)
                    call date_and_time(date, time, zone, values)
!                   write(*,FMT='(A10,A,I4,A,I4)') time, ' - Master ', myrank, ' - Read observation ', n_procs_done
                    
!                   write(24,FMT='(A10,A,I4,A,I4)') time, ' Master ', myrank, ' - Read observation ', n_procs_done
                    if (status_obs == 1) then
                        call send_observation(observation, package_size_obs, (n_procs_done - starting_pixel + 1), n_procs_done)
                        call date_and_time(date, time, zone, values)
                        write(*,FMT='(A10,A,I4,A,I5,A,I4)') time, ' * Master - Sent observation ', n_procs_done, ' to slave ', (n_procs_done - starting_pixel + 1)
                        
!                       write(24,FMT='(A10,A,I4,A,I5,A,I4)') time, ' Master ', myrank, ' - Sent observation ', n_procs_done, ' to Slave ', &
!                           (n_procs_done - starting_pixel + 1)
                        n_procs_done = n_procs_done + 1
                    endif                   
                enddo               

! If the number of observations is smaller than the number of slaves, kill the remaining slaves             
                if ((n_procs_done - starting_pixel + 1) < nprocs) then
                    do slave = (n_procs_done - starting_pixel + 1), nprocs-1
                        call kill_slave(slave)
                        write(*,FMT='(A,I4)') 'Too many slaves. Kill signal sent to slave ', slave              
!                       write(24,FMT='(A,I4)') 'Too many slaves. Kill signal sent to Slave ', slave
                        slave_active(slave) = 0
                    enddo
                endif

!               call Flush(24)
                
            endif
            
! Do while we have something new to do
            do while (status_obs == 1 .or. sum(slave_active) /= 0)
            
! If master, receive final model, read new observation and send it to the slave that finished the work
                if (myrank == 0) then               

                    call receive_model(params, errorparams, observation, inversion, package_size_model, slave, index_obs)
!                   call date_and_time(date, time, zone, values)
!                   write(*,FMT='(A10,A,I4,A,I4)') time, ' - Master ', myrank, ' - Received result from Slave ', slave
                    
!                   write(24,FMT='(A10,A,I4,A,I4)') time, ' Master ', myrank, ' - Received result from Slave ', slave

                    call write_results(fixed, observation, inversion, params, errorparams, index_obs)
                    call date_and_time(date, time, zone, values)
                    write(*,FMT='(A10,A,I5,A,I4)') time, ' * Master - Received and saved result ', index_obs, ' from slave ', slave
                    
                    
                    call read_observation(input_observed_profiles, fixed, observation, params, n_procs_done, status_obs)

                    if (status_obs == 1) then
                        call send_observation(observation, package_size_obs, slave, n_procs_done)
                        call date_and_time(date, time, zone, values)
                        write(*,FMT='(A10,A,I5,A,I4)') time, ' * Master - Sent observation ', n_procs_done, ' to slave ', slave
                        
                    else
! No more observations, so send a kill order to the slave
                        if (slave_active(slave) == 1) then
                            call kill_slave(slave)
                            call date_and_time(date, time, zone, values)
                            write(*,FMT='(A10,A,I4)') time, ' * Master - Kill signal sent to slave ', slave
                            
                            slave_active(slave) = 0
                        endif
                    endif
                    

                    n_procs_done = n_procs_done + 1

!                   call Flush(24)
                    
                endif

! If slave, receive observations, do the work and send the new inferred model to the master
                if (myrank /= 0) then

                    call receive_observation(observation, package_size_obs,kill,index_obs)

! If we received a new observation, do the inversion
                    if (kill == 0) then
!                       write(*,FMT='(A,I4,A)') 'Slave  ', myrank, ' - Received observation'
!                       write(*,FMT='(A,I4,A)') 'Slave  ', myrank, ' - Doing inversion'
                        call set_boundary_condition(fixed, inversion)
                        
                        if (working_mode >= 1) then
                            call doinversion(params, errorparams, fixed, observation, inversion, myrank, error)
                        endif

                        if (working_mode == 0) then
                            call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
                        endif

                        if (error == 0) then
                            call send_model(params, errorparams, observation, inversion, package_size_model, myrank, index_obs)
                        else
!                           write(*,FMT='(A,I4,A)') 'Slave  ', myrank, ' - Model not converged'
                        endif
!                       write(*,FMT='(A,I4,A)') 'Slave  ', myrank, ' - Sent back result'
                    else

! If not, quit
                        status_obs = 0
                        
                    endif                   
                endif
                
            enddo

        else

! If only one node is available, loop over all observations and carry
! out the inversion serially
            status_obs = 1
            do while (status_obs == 1)
                call read_observation(input_observed_profiles, fixed, observation, params, n_procs_done, status_obs)                
                call set_boundary_condition(fixed, inversion)
                call doinversion(params, errorparams, fixed, observation, inversion, myrank, error)
            enddo
            
        endif

! Close files
        if (myrank == 0) then                           
#if defined(mpi)
            call check( nf90_close(observation%obs_id) )
            call check( nf90_close(fixed%syn_id) )
            call check( nf90_close(fixed%par_id) )
            call check( nf90_close(fixed%error_id) )
#endif
        endif
    
    endif
    
    if (myrank == 0) then
        open(unit=12,file='done.info',action='write',status='replace')
        write(12,*) 'Done.'
        close(12)
    endif

#if defined(mpi)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)    
    Call MPI_FINALIZE(mpi_status)
#endif
        
end program multiterm