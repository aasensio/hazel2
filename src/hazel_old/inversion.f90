module inversion_mod
use vars, only : variable_parameters, fixed_parameters, type_observation, type_inversion
use marquardt
use maths, only : secantConfidenceLevel
implicit none
contains

!------------------------------------------------------------
! Carry out the inversion
! If myrank /= 0, then we are using a master/slave strategy and we output
! only relevant information
!------------------------------------------------------------
    subroutine doinversion(params, errorparams, fixed, observation, inversion, myrank, error)
    type(variable_parameters) :: params, errorparams
    type(fixed_parameters) :: fixed
    type(type_observation) :: observation
    type(type_inversion) :: inversion
    integer :: myrank, error

    integer :: iter, nworst_chisq, loop_cycle, successful
    real(kind=8) :: chisq_relative_change, params_relative_change, nu
    logical :: correct

        error = 0
        
! Compute confidence level depending on the number of degrees of freedom
! We use an approximation to the number of degrees of freedom using the total number of parameters
! because this is dominated by the number of observations
        nu = 4.d0 * observation%n - params%n_total
        chi2Level = secantConfidenceLevel(nu, erf(1.d0/sqrt(2.d0)))
                
! If the first cycle is LM, carry out a first synthesis
! with the original values of the parameters
        if (inversion%algorithm(1) == 1) then
            inversion%loop_cycle = 1        
            call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
            if (error == 1) return
            inversion%chisq = compute_chisq(observation,inversion)
        endif
        
        errorparams = params
        
! Loop over the number of cycles
        do loop_cycle = 1, inversion%n_cycles
            
            inversion%loop_cycle = loop_cycle
                
! Select the parameters to invert in this cycle
            params%inverted = inversion%cycles(:,loop_cycle)
                            
            if (myrank == 0) then
                write(*,FMT='(A)') '*******************************'
                write(*,FMT='(A,I2,A1,I2)') 'Starting cycle ', inversion%loop_cycle, '/', inversion%n_cycles
                write(*,FMT='(A,I2,A)') ' Inverting ', sum(params%inverted),' parameters'
                write(*,FMT='(A)') '*******************************'
                
                write(*,FMT='(A,4(2X,F8.3))') 'Stokes parameters weights : ', inversion%stokes_weights(0:3,loop_cycle)
            endif
    
! LEVENBERG-MARQUARDT
            if (inversion%algorithm(loop_cycle) == 1) then
    
! If the weights of the previous step were different, recalculate the value of chi^2
                if (loop_cycle > 1) then
                    if (inversion%stokes_weights(0,loop_cycle-1) /= inversion%stokes_weights(0,loop_cycle) .or. &
                        inversion%stokes_weights(1,loop_cycle-1) /= inversion%stokes_weights(1,loop_cycle) .or. &
                        inversion%stokes_weights(2,loop_cycle-1) /= inversion%stokes_weights(2,loop_cycle) .or. &
                        inversion%stokes_weights(3,loop_cycle-1) /= inversion%stokes_weights(3,loop_cycle)) then

                        call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)                      
                        if (error == 1) return
                        
                        inversion%chisq = compute_chisq(observation,inversion)
                    endif
                endif

                if (myrank == 0) then
                    print *, 'LEVENBERG-MARQUARDT MODE'
                endif
    
                inversion%lambda = 0.001d0
                nworst_chisq = 0
                inversion%min_chisq = inversion%chisq
                inversion%chisq = 1.d10
                inversion%chisq_old = inversion%chisq           
                chisq_relative_change = 1.d10
                iter = 1
                successful = 1

                call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
                if (error == 1) return
    
! Main inversion loop
                do while (iter < inversion%iter_max .and. nworst_chisq < 10 .and. abs(chisq_relative_change) > 1.d-4 )
    
                    if (myrank == 0) then
                        write(*,FMT='(A,I3,A1,I3,A,F11.6)') 'Iteration ', iter, '/', inversion%iter_max, ' -- Lambda : ', inversion%lambda
                    endif

! Only recalculate derivatives if the step has been successful (we located a point with a smaller chi^2)
                    if (successful /= 0) call compute_dydx(params,fixed,inversion,observation,error)

                    if (error == 1) return

                    call compute_trial_params(params,fixed,inversion,observation,trial)
    
                    call check_boundaries(params,trial,correct,fixed)

                    if (myrank == 0) then
                        call print_parameters(params,'  -Old parameters : ',.TRUE.)
                        call print_parameters(trial,'  -New parameters : ',.FALSE.)
                    endif

!                   write(*,FMT='(10X,A,I4,A,F18.8)') 'LM - chi^2(', myrank, ') : ', inversion%chisq
    
                    if (correct) then
                        call do_synthesis(trial, fixed, observation, inversion%stokes_unperturbed, error)
                        if (error == 1) return
                        
                        inversion%chisq_old = inversion%chisq
                        inversion%chisq = compute_chisq(observation,inversion)
                        
                        call compute_uncertainty(trial,fixed,inversion,observation,errorparams)
    
! Verify if the chisq is smaller or larger than the previous fit
! Better model
                        if (inversion%chisq < inversion%min_chisq) then
                            successful = 1                          

! Change lambda
                            if (inversion%lambda >= 1.d4) then
                                inversion%lambda = inversion%lambda / 100.d0
                            else if (inversion%lambda >= 1.d-4 .and. inversion%lambda < 1.d4) then
                                inversion%lambda = inversion%lambda / 10.d0
                            else if (inversion%lambda < 1.d-4) then
                                inversion%lambda = inversion%lambda / 5.d0
                            endif
                            
                            if (inversion%lambda < 1.d-6) inversion%lambda = 1.d-6
                            !params_relative_change = compute_params_relative_change(params,trial)
                            params = trial
                            chisq_relative_change = (inversion%chisq-inversion%chisq_old) / inversion%chisq * 100.d0
                            if (myrank == 0) then       
                                print *, 'New chi^2 : ', inversion%chisq                            
                                print *, 'chi^2 relative change [%] : ', chisq_relative_change
                                print *, 'Relative change in the parameters [%] : ', params_relative_change
                                print *, 'Successful step'
                            endif
                            nworst_chisq = 0
!                           call write_final_profiles(temporal_file,observation,inversion)
                            inversion%min_chisq = inversion%chisq
                            !successful = successful + 1
                        else
! Worse model
                            successful = 0

                            if (inversion%lambda >= 1.d4) then
                                inversion%lambda = inversion%lambda * 100.d0
                            else if (inversion%lambda >= 1.d-4 .and. inversion%lambda < 1.d4) then
                                inversion%lambda = inversion%lambda * 10.d0
                            else if (inversion%lambda < 1.d-4) then
                                inversion%lambda = inversion%lambda * 5.d0
                            endif

                            nworst_chisq = nworst_chisq + 1
                            if (myrank == 0) then
                                print *, 'Larger chi^2 : ', inversion%chisq, ' --- Increasing Lambda'                           
                                write(*,FMT='(A,I1,A1,I1)') '   Step : ', nworst_chisq, '/', 5
                            endif
    
                        endif           
                    else
                        if (myrank == 0) then
                            print *, 'Unphysical values. Trying again...'
                        endif
                        inversion%lambda = inversion%lambda * 10.d0
                    endif

                    if (myrank == 0) then
                        print *, '---------------------'
                        print *
                    endif
    
                    iter = iter + 1
    
                enddo  ! Iterations

                if (myrank == 0) then
                    print *, 'Number of calls to forward modeling routines : ', fixed%total_forward_modeling
                endif
    
                call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
                if (error == 1) return
                inversion%chisq = compute_chisq(observation,inversion)
                
                if (myrank == 0) then

! Write the final profiles
                    call write_final_profiles(output_inverted_profiles,observation,inversion)

! Write the final parameters in a file so that it can be used for restarting the inversion code
                    call write_experiment(params, fixed, errorparams)
                
                    call print_parameters(params,'-Final Parameters : ',.TRUE.)
                    print *, 'Final chi^2 : ', inversion%chisq
                endif
                
            endif
                
                    
!*********************************
!** INVERSION MODE WITH DIRECT
!*********************************
            if (inversion%algorithm(loop_cycle) == 2 .or. inversion%algorithm(loop_cycle) == 3) then
                if (myrank == 0) then
                    print *, 'DIRECT MODE'
                endif
                    
                trial = params
                call invert_with_direct(trial,fixed, myrank, inversion%algorithm(loop_cycle), error)

                if (error == 1) return
                params = trial
                
                call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
                
                if (error == 1) return
                inversion%chisq = compute_chisq(observation,inversion)              

                if (myrank == 0) then
! Write the final profiles
                    call write_final_profiles(output_inverted_profiles,observation,inversion)
                
! Write the final parameters in a file so that it can be used for restarting the inversion code
                    call write_experiment(params, fixed, errorparams)
                    
                    call print_parameters(params,'-Final Parameters : ',.TRUE.)
                    print *, 'Final chi^2 : ', inversion%chisq
                endif
                            
            endif
        
! !*********************************
! !** INVERSION MODE WITH PIKAIA
! !*********************************
!           if (inversion%algorithm(loop_cycle) == 3) then
!               if (myrank == 0) then
!                   print *, 'PIKAIA MODE'
!               endif
!   
!               trial = params
!               call invert_with_pikaia(trial)
!               params = trial
!               
!               call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)
! 
!               if (error == 1) return
!               
!               inversion%chisq = compute_chisq(observation,inversion)
! 
!               if (myrank == 0) then
! ! Write the final profiles
!                   call write_final_profiles(output_inverted_profiles,observation,inversion)
!               
! ! Write the final parameters in a file so that it can be used for restarting the inversion code
!                   call write_experiment(params, fixed)
!               
!                   call print_parameters(params,'-Final Parameters : ',.TRUE.)
!                   print *, 'Final chi^2 : ', inversion%chisq
!               endif
!           endif

            if (myrank == 0) then
                print *
                write(*,FMT='(A)') '*******************************'
                write(*,FMT='(A,I2,A1,I2)') 'End of cycle ', inversion%loop_cycle, '/', inversion%n_cycles
                write(*,FMT='(A)') '*******************************'
                print *            
            endif
                
        enddo  ! Cycles

    end subroutine doinversion

end module inversion_mod
