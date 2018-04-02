module mpi_routines_fake
use vars
implicit none

contains

!------------------------------------------------------------
! Fake MPI barrier
!------------------------------------------------------------
    subroutine MPI_Barrier(world, ierr)
    integer :: world, ierr
        
    end subroutine MPI_Barrier

!------------------------------------------------------------
! Fake MPI barrier
!------------------------------------------------------------
    subroutine MPI_Finalize(ierr)
    integer :: ierr
        
    end subroutine MPI_Finalize

!------------------------------------------------------------
! Broadcast Allen's data
!------------------------------------------------------------
    subroutine bcast_allen
        
    end subroutine bcast_allen

!------------------------------------------------------------
! Broadcast main configuration
!------------------------------------------------------------
    subroutine bcast_main_conf
        
    end subroutine bcast_main_conf

!------------------------------------------------------------
! Broadcast experiment
!------------------------------------------------------------
    subroutine bcast_experiment_conf
        
    end subroutine bcast_experiment_conf


!------------------------------------------------------------
! Broadcast atom model
!------------------------------------------------------------
    subroutine bcast_atom_model(myrank)
    integer :: myrank
       
    end subroutine bcast_atom_model

!------------------------------------------------------------
! Broadcast main configuration
!------------------------------------------------------------
    subroutine bcast_parameters_to_invert(myrank)
    integer :: myrank
        
    end subroutine bcast_parameters_to_invert

!------------------------------------------------------------
! Broadcast information about sizes of observations
!------------------------------------------------------------
    subroutine bcast_observation_info
        
    end subroutine bcast_observation_info

!------------------------------------------------------------
! Broadcast DIRECT configuration
!------------------------------------------------------------
    subroutine bcast_direct_info(myrank)
    integer :: myrank

    end subroutine bcast_direct_info
!------------------------------------------------------------
! Send a new observation to a slave
!------------------------------------------------------------
    subroutine send_observation(observation, packagesize, slave, index_obs)
    type(type_observation) :: observation
    integer :: slave, packagesize, pos, index_obs
            
    end subroutine send_observation

!------------------------------------------------------------
! Receive a new observation from the master
!------------------------------------------------------------
    subroutine receive_observation(observation, packagesize, stop_flag, index_obs)
    type(type_observation) :: observation
    integer :: packagesize, pos, stop_flag, index_obs
        
    end subroutine receive_observation

!------------------------------------------------------------
! Send a converged model to the master
!------------------------------------------------------------
    subroutine send_model(params, error, in_observation, in_inversion, packagesize, myrank, index_obs)
    type(variable_parameters) :: params, error
    type(type_observation) :: in_observation
    type(type_inversion) :: in_inversion
    integer :: packagesize, pos, myrank, index_obs
    
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
            
    end subroutine receive_model
    

!------------------------------------------------------------
! Send the kill signal to all slaves
!------------------------------------------------------------
    subroutine kill_slave(slave)
    integer :: slave
        
    end subroutine kill_slave

end module mpi_routines_fake