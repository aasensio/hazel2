module marquardt
use vars
use SEE
use rt_coef
use synth
use svd
use maths
use io
implicit none
contains

!*************************************************************
!*************************************************************
! MARQUARDT METHOD
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Calculate the chi^2 function
!------------------------------------------------------------
	function compute_chisq(in_observation,in_inversion)
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	real(kind=8) :: compute_chisq, weight
	integer :: nf, i, j
		
		compute_chisq = 0.d0
		
! Compute the chisq using the appropriate weights for each cycle
		do i = 0, 3
			weight = in_inversion%stokes_weights(i,in_inversion%loop_cycle)
			compute_chisq = compute_chisq + weight * &
				sum((in_observation%stokes(i,:)-in_inversion%stokes_unperturbed(i,:))**2/in_observation%sigma(i,:)**2) / &
				(4.d0*in_observation%n)			
		enddo
		
	end function compute_chisq
	
!------------------------------------------------------------
! Calculate the maximum relative change in the parameters
!------------------------------------------------------------
	function compute_params_relative_change(params,trial)
	real(kind=8) :: compute_params_relative_change
	type(variable_parameters) :: params, trial
	real(kind=8), allocatable :: rel_change(:)
	integer :: i, j
	
		allocate(rel_change(params%n_total))

		if (params%bgauss /= 0.d0) then
			rel_change(1) = abs((params%bgauss-trial%bgauss) / params%bgauss)
		else
			rel_change(1) = 0.d0
		endif
		
		if (params%thetabd /= 0.d0) then
			rel_change(2) = abs((params%thetabd-trial%thetabd) / params%thetabd)
		else
			rel_change(2) = 0.d0
		endif
		
		if (params%chibd /= 0.d0) then
			rel_change(3) = abs((params%chibd-trial%chibd) / params%chibd)
		else
			rel_change(3) = 0.d0
		endif
		
		if (params%vdopp /= 0.d0) then
			rel_change(4) = abs((params%vdopp-trial%vdopp) / params%vdopp)
		else
			rel_change(4) = 0.d0
		endif
		
		if (params%dtau /= 0.d0) then
			rel_change(5) = abs((params%dtau-trial%dtau) / params%dtau)
		else
			rel_change(5) = 0.d0
		endif
		
		if (params%delta_collision /= 0.d0) then
			rel_change(6) = abs((params%delta_collision-trial%delta_collision) / params%delta_collision)
		else
			rel_change(6) = 0.d0
		endif
		
		if (params%vmacro /= 0.d0) then
			rel_change(7) = abs((params%vmacro-trial%vmacro) / params%vmacro)
		else
			rel_change(7) = 0.d0
		endif
		
		if (params%damping /= 0.d0) then
			rel_change(8) = abs((params%damping-trial%damping) / params%damping)
		else
			rel_change(8) = 0.d0
		endif
		
		if (params%beta /= 0.d0) then
			rel_change(9) = abs((params%beta-trial%beta) / params%beta)
		else
			rel_change(9) = 0.d0
		endif
		
		if (params%height /= 0.d0) then
			rel_change(10) = abs((params%height-trial%height) / params%height)
		else
			rel_change(10) = 0.d0
		endif

		if (params%dtau2 /= 0.d0) then
			rel_change(11) = abs((params%dtau2-trial%dtau2) / params%dtau2)
		else
			rel_change(11) = 0.d0
		endif

		if (params%vmacro2 /= 0.d0) then
			rel_change(12) = abs((params%vmacro2-trial%vmacro2) / params%vmacro2)
		else
			rel_change(12) = 0.d0
		endif

		if (params%bgauss2 /= 0.d0) then
			rel_change(13) = abs((params%bgauss2-trial%bgauss2) / params%bgauss2)
		else
			rel_change(13) = 0.d0
		endif
		
		if (params%thetabd2 /= 0.d0) then
			rel_change(14) = abs((params%thetabd2-trial%thetabd2) / params%thetabd2)
		else
			rel_change(14) = 0.d0
		endif
		
		if (params%chibd2 /= 0.d0) then
			rel_change(15) = abs((params%chibd2-trial%chibd2) / params%chibd2)
		else
			rel_change(15) = 0.d0
		endif

		if (params%vdopp2 /= 0.d0) then
			rel_change(16) = abs((params%vdopp2-trial%vdopp2) / params%vdopp2)
		else
			rel_change(16) = 0.d0
		endif

		if (params%ff /= 0.d0) then
			rel_change(17) = abs((params%ff-trial%ff) / params%ff)
		else
			rel_change(17) = 0.d0
		endif

		if (params%beta2 /= 0.d0) then
			rel_change(18) = abs((params%beta-trial%beta) / params%beta)
		else
			rel_change(18) = 0.d0
		endif
		
		compute_params_relative_change = maxval(rel_change)*100.d0		

		deallocate(rel_change)
		
		
	end function compute_params_relative_change
	
	
!------------------------------------------------------------
! Check if the values are inside physical boundaries
! If not inside the boundaries, use the previous value
!------------------------------------------------------------
	subroutine check_boundaries(original,trial,correct,fixed)
	type(variable_parameters) :: original, trial
	type(fixed_parameters) :: fixed
	logical :: correct
	real(kind=8) :: upper(18), lower(18)
	integer :: i
		
! Verify if the parameters are inside the boundaries. If any of them is outside, use the previous value
! of this parameter.
	
		if (trial%bgauss < fixed%lower_direct(1) .or. trial%bgauss > fixed%upper_direct(1)) then
			trial%bgauss = original%bgauss
		endif
		if (trial%thetabd < fixed%lower_direct(2) .or. trial%thetabd > fixed%upper_direct(2)) then
			trial%thetabd = original%thetabd
		endif
		if (trial%chibd < fixed%lower_direct(3) .or. trial%chibd > fixed%upper_direct(3)) then
			trial%chibd = original%chibd
		endif
		if (trial%vdopp < fixed%lower_direct(4) .or. trial%vdopp > fixed%upper_direct(4)) then
			trial%vdopp = original%vdopp
		endif
		if (trial%dtau < fixed%lower_direct(5) .or. trial%dtau > fixed%upper_direct(5)) then
			trial%dtau = original%dtau
		endif
		if (original%nslabs == 2) then
			if (trial%dtau2 < fixed%lower_direct(11) .or. trial%dtau2 > fixed%upper_direct(11)) then
				trial%dtau2 = original%dtau2
			endif
		endif
		if (trial%delta_collision < fixed%lower_direct(6) .or. trial%delta_collision > fixed%upper_direct(6)) then
			trial%delta_collision = original%delta_collision
		endif		
		if (trial%vmacro < fixed%lower_direct(7) .or. trial%vmacro > fixed%upper_direct(7)) then
			trial%vmacro = original%vmacro
		endif
		if (original%nslabs == 2) then
			if (trial%vmacro2 < fixed%lower_direct(12) .or. trial%vmacro2 > fixed%upper_direct(12)) then
				trial%vmacro2 = original%vmacro2
			endif
		endif
		if (trial%damping < fixed%lower_direct(8) .or. trial%damping > fixed%upper_direct(8)) then
			trial%damping = original%damping
		endif
		if (trial%beta < fixed%lower_direct(9) .or. trial%beta > fixed%upper_direct(9)) then
			trial%beta = original%beta
		endif
		if (trial%height < fixed%lower_direct(10) .or. trial%height > fixed%upper_direct(10)) then
			trial%height = original%height
		endif

! If two components with two different fields
		if (original%nslabs == 3) then
			if (trial%bgauss2 < fixed%lower_direct(13) .or. trial%bgauss2 > fixed%upper_direct(13)) then
				trial%bgauss2 = original%bgauss2
			endif
			if (trial%thetabd2 < fixed%lower_direct(14) .or. trial%thetabd2 > fixed%upper_direct(14)) then
				trial%thetabd2 = original%thetabd2
			endif
			if (trial%chibd2 < fixed%lower_direct(15) .or. trial%chibd2 > fixed%upper_direct(15)) then
				trial%chibd2 = original%chibd2
			endif
			if (trial%vdopp2 < fixed%lower_direct(16) .or. trial%vdopp2 > fixed%upper_direct(16)) then
				trial%vdopp2 = original%vdopp2
			endif
		endif

! Two slabs in the same pixel
		if (original%nslabs == -2) then
			if (trial%ff < fixed%lower_direct(17) .or. trial%ff > fixed%upper_direct(17)) then
				trial%ff = original%ff
			endif
		endif

		if (original%nslabs == 3 .or. original%nslabs == -2 .or. original%nslabs == -2) then
			if (trial%beta2 < fixed%lower_direct(18) .or. trial%beta2 > fixed%upper_direct(18)) then
				trial%beta2 = original%beta2
			endif
		endif

		
		correct = .TRUE.
		
	end subroutine check_boundaries
	

!------------------------------------------------------------
! Perturb a given parameter
!------------------------------------------------------------
	subroutine perturb_parameter(in_params,out_params,which,perturbation,delta)
	type(variable_parameters) :: in_params, out_params
	integer :: which
	real(kind=8) :: perturbation, delta

		select case(which)
			case(1) 
				if (in_params%bgauss > 1.d-2) then
					out_params%bgauss = in_params%bgauss * (1.d0+perturbation)
				else
					out_params%bgauss = in_params%bgauss + perturbation
				endif
				delta = out_params%bgauss - in_params%bgauss
				
			case(2)
				if (abs(in_params%thetabd) > 1.d-2) then
					out_params%thetabd = in_params%thetabd * (1.d0+perturbation)
				else
					out_params%thetabd = in_params%thetabd + perturbation
				endif
				delta = out_params%thetabd - in_params%thetabd
				
			case(3) 
				if (abs(in_params%chibd) > 1.d-2) then				
					out_params%chibd = in_params%chibd * (1.d0+perturbation)
				else
					out_params%chibd = in_params%chibd + perturbation
				endif
				delta = out_params%chibd - in_params%chibd
				
			case(4) 
				if (in_params%vdopp > 1.d-2) then
					out_params%vdopp = in_params%vdopp * (1.d0+perturbation)
				else
					out_params%vdopp = in_params%vdopp + perturbation
				endif
				delta = out_params%vdopp - in_params%vdopp
				
			case(5) 
				if (in_params%dtau > 1.d-2) then
					out_params%dtau = in_params%dtau * (1.d0+perturbation)
				else
					out_params%dtau = in_params%dtau + perturbation
				endif
				delta = out_params%dtau - in_params%dtau
				
			case(6) 
				if (in_params%delta_collision > 1.d-2) then
					out_params%delta_collision = in_params%delta_collision * (1.d0+perturbation)
				else
					out_params%delta_collision = in_params%delta_collision + perturbation
				endif
				delta = out_params%delta_collision - in_params%delta_collision
				
			case(7)
				if (abs(in_params%vmacro) > 1.d-2) then				
					out_params%vmacro = in_params%vmacro * (1.d0+perturbation)
				else
					out_params%vmacro = in_params%vmacro + perturbation
				endif
				delta = out_params%vmacro - in_params%vmacro
				
			case(8)
				if (in_params%damping > 1.d-2) then
					out_params%damping = in_params%damping * (1.d0+perturbation)
				else
					out_params%damping = in_params%damping + perturbation
				endif
				delta = out_params%damping - in_params%damping
				
			case(9)
				if (abs(in_params%beta) > 1.d-2) then
					out_params%beta = in_params%beta * (1.d0+perturbation)
				else
					out_params%beta = in_params%beta + perturbation
				endif
				delta = out_params%beta - in_params%beta
				
			case(10)
				if (in_params%height > 1.d-2) then
					out_params%height = in_params%height * (1.d0+perturbation)
				else
					out_params%height = in_params%height + perturbation
				endif
				delta = out_params%height - in_params%height

			case(11) 
				if (in_params%dtau2 > 1.d-2) then
					out_params%dtau2 = in_params%dtau2 * (1.d0+perturbation)
				else
					out_params%dtau2 = in_params%dtau2 + perturbation
				endif
				delta = out_params%dtau2 - in_params%dtau2

			case(12)
				if (abs(in_params%vmacro2) > 1.d-2) then				
					out_params%vmacro2 = in_params%vmacro2 * (1.d0+perturbation)
				else
					out_params%vmacro2 = in_params%vmacro2 + perturbation
				endif
				delta = out_params%vmacro2 - in_params%vmacro2

			case(13) 
				if (in_params%bgauss2 > 1.d-2) then
					out_params%bgauss2 = in_params%bgauss2 * (1.d0+perturbation)
				else
					out_params%bgauss2 = in_params%bgauss2 + perturbation
				endif
				delta = out_params%bgauss2 - in_params%bgauss2
				
			case(14)
				if (abs(in_params%thetabd2) > 1.d-2) then
					out_params%thetabd2 = in_params%thetabd2 * (1.d0+perturbation)
				else
					out_params%thetabd2 = in_params%thetabd2 + perturbation
				endif
				delta = out_params%thetabd2 - in_params%thetabd2
				
			case(15) 
				if (abs(in_params%chibd2) > 1.d-2) then				
					out_params%chibd2 = in_params%chibd2 * (1.d0+perturbation)
				else
					out_params%chibd2 = in_params%chibd2 + perturbation
				endif
				delta = out_params%chibd2 - in_params%chibd2

			case(16) 
				if (in_params%vdopp2 > 1.d-2) then
					out_params%vdopp2 = in_params%vdopp2 * (1.d0+perturbation)
				else
					out_params%vdopp2 = in_params%vdopp2 + perturbation
				endif
				delta = out_params%vdopp2 - in_params%vdopp2

			case(17)
				if (in_params%ff > 1.d-2) then
					out_params%ff = in_params%ff * (1.d0+perturbation)
				else
					out_params%ff = in_params%ff + perturbation
				endif
				delta = out_params%ff - in_params%ff

			case(18)
				if (abs(in_params%beta2) > 1.d-2) then
					out_params%beta2 = in_params%beta2 * (1.d0+perturbation)
				else
					out_params%beta2 = in_params%beta2 + perturbation
				endif
				delta = out_params%beta2 - in_params%beta2
			
		end select
		
	end subroutine perturb_parameter
	
!------------------------------------------------------------
! Perturb a given parameter
!------------------------------------------------------------
	subroutine add_to_parameter(in_params,out_params,which,new_value)
	type(variable_parameters) :: in_params, out_params
	integer :: which
	real(kind=8) :: new_value, rel_change

! Limit the maximum relative change to 50%
		if (new_value > 0.5d0) then
! 			write(*,FMT='(A,A)') 'Clipping ', trim(adjustl(parameters_name(which)))
			new_value = 0.5d0
		endif
		select case(which)
			case(1)				
				out_params%bgauss = in_params%bgauss + new_value
			case(2) 
				out_params%thetabd = in_params%thetabd + new_value
			case(3) 
				out_params%chibd = in_params%chibd + new_value
			case(4) 
				out_params%vdopp = in_params%vdopp + new_value
			case(5) 
				out_params%dtau = in_params%dtau + new_value
			case(6) 
				out_params%delta_collision = in_params%delta_collision + new_value
			case(7) 
				out_params%vmacro = in_params%vmacro + new_value	
			case(8) 
				out_params%damping = in_params%damping + new_value
			case(9) 
				out_params%beta = in_params%beta + new_value
			case(10) 
				out_params%height = in_params%height + new_value
			case(11) 
				out_params%dtau2 = in_params%dtau2 + new_value
			case(12) 
				out_params%vmacro2 = in_params%vmacro2 + new_value
			case(13)				
				out_params%bgauss2 = in_params%bgauss2 + new_value
			case(14) 
				out_params%thetabd2 = in_params%thetabd2 + new_value
			case(15) 
				out_params%chibd2 = in_params%chibd2 + new_value
			case(16) 
				out_params%vdopp2 = in_params%vdopp2 + new_value
			case(17)
				out_params%ff = in_params%ff + new_value
			case(18) 
				out_params%beta2 = in_params%beta2 + new_value
		end select
		
	end subroutine add_to_parameter	
	
!------------------------------------------------------------
! Calculate the derivatives of the Stokes profiles
!------------------------------------------------------------
	subroutine compute_dydx(in_params,in_fixed,in_inversion,in_observation,error)
	type(variable_parameters) :: in_params, params_scaled, params_scaled_modified, params_modified
	type(fixed_parameters) :: in_fixed
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	integer :: i, error
	real(kind=8) :: delta
								
! Do the synthesis without perturbing the magnetic field		
! 		call do_synthesis(in_params, in_fixed, in_observation, in_inversion%stokes_unperturbed)
		
		params_scaled = in_params		
		params_modified = in_params

		call compress_parameters(in_params, params_scaled)		
		params_scaled_modified = params_scaled

		do i = 1, in_params%n_total
		
			if (in_params%inverted(i) == 1) then
!  				 write(*,FMT='(A,I2,A,A)') 'Perturbing parameter ', i, '  --- ', trim(adjustl(parameters_name(i)))
		
! Perturb the given parameter and do another synthesis
				call perturb_parameter(params_scaled,params_scaled_modified,i,0.0001d0,delta)
				
				call expand_parameters(params_scaled_modified, params_modified)

				call do_synthesis(params_modified, in_fixed, in_observation, in_inversion%stokes_perturbed, error)

				if (error == 1) return
		
				in_inversion%dydx(:,:,i) = (in_inversion%stokes_perturbed - in_inversion%stokes_unperturbed) / delta

			endif
		enddo
	
	end subroutine compute_dydx
	
!------------------------------------------------------------
! Calculate the derivatives of the Stokes profiles
!------------------------------------------------------------
	subroutine compute_trial_params(in_params,in_fixed,in_inversion,in_observation,in_trial)
	type(variable_parameters) :: in_params, in_trial, in_scaled, in_temp
	type(fixed_parameters) :: in_fixed
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	real(kind=8), allocatable :: alpha(:,:), beta(:), w(:), v(:,:), x(:)
	real(kind=8) :: obs, syn, sig, weight, wmax
	integer :: i, j, k, l, np
			
		in_trial = in_params		
		call compress_parameters(in_params,in_scaled)
		in_temp = in_scaled
				
		np = in_params%n_total
		
		allocate(alpha(np,np))
		allocate(beta(np))
		allocate(v(np,np))
		allocate(w(np))
		allocate(x(np))
		
		alpha = 0.d0
		beta = 0.d0
						
		do i = 1, in_observation%n
			do j = 0, 3
				weight = in_inversion%stokes_weights(j,in_inversion%loop_cycle)
				do k = 1, np
					
					obs = in_observation%stokes(j,i)
					syn = in_inversion%stokes_unperturbed(j,i)
					sig = in_observation%sigma(j,i)
					
					if (in_params%inverted(k) == 1) then
						beta(k) = beta(k) - 2.d0 * weight * (obs-syn) / sig**2 * in_inversion%dydx(j,i,k) / (4.d0*in_observation%n)
						do l = 1, np
							if (in_params%inverted(l) == 1) then
								alpha(k,l) = alpha(k,l) + 2.d0 * weight * &
									in_inversion%dydx(j,i,k) * in_inversion%dydx(j,i,l) / sig**2 / (4.d0*in_observation%n)
							endif
						enddo
					endif
				enddo
			enddo
		enddo
		
		beta = -0.5d0 * beta
		alpha = 0.5d0 * alpha
		
! Apply the Lambda parameter of the Levenverb-Marquardt		
		do i = 1, np
			alpha(i,i) = (1.d0+inversion%lambda) * alpha(i,i)
		enddo

		np = in_params%n_total
		call svdcmp(alpha,np,np,np,np,w,v)

! Regularize the linear system by setting all singular values below a certain
! threshold equal to zero
		wmax = maxval(w)
		do i = 1, np
			if (w(i) < wmax * 1.d-6) then
				w(i) = 0.d0
			endif
		enddo
		
		call svbksb(alpha,w,v,np,np,np,np,beta,x)				

		do i = 1, np
			if (in_params%inverted(i) == 1) then								
				call add_to_parameter(in_scaled,in_temp,i,x(i)/1.d0)
			endif				
		enddo
				
		call expand_parameters(in_temp,in_trial)
		
		deallocate(alpha)
		deallocate(beta)
		deallocate(v)
		deallocate(w)
		deallocate(x)
	
	end subroutine compute_trial_params
	
!------------------------------------------------------------
! Calculate the derivatives of the Stokes profiles
!------------------------------------------------------------
	subroutine compute_uncertainty(in_params,in_fixed,in_inversion,in_observation,error)
	type(variable_parameters) :: in_params, in_scaled, in_temp, error
	type(fixed_parameters) :: in_fixed
	type(type_inversion) :: in_inversion
	type(type_observation) :: in_observation
	real(kind=8), allocatable :: alpha(:,:), beta(:,:), w(:), v(:,:), x(:)
	real(kind=8) :: obs, syn, sig, weight, wmax, value
	integer :: i, j, k, l, np
					
		call compress_parameters(in_params,in_scaled)
		in_temp = in_scaled
				
		np = in_params%n_total
				
		allocate(alpha(np,np))		
		allocate(beta(np,np))
		allocate(v(np,np))
		allocate(w(np))
		allocate(x(np))
		
		alpha = 0.d0
		beta = 0.d0
								
		do i = 1, in_fixed%no
			do j = 0, 3				
				do k = 1, np
					
					obs = in_observation%stokes(j,i)
					syn = in_inversion%stokes_unperturbed(j,i)
					sig = in_observation%sigma(j,i)
															
					if (in_params%inverted(k) == 1) then
						do l = 1, np
							if (in_params%inverted(l) == 1) then
								alpha(k,l) = alpha(k,l) + 2.d0 * in_inversion%dydx(j,i,k) * in_inversion%dydx(j,i,l) / sig**2
							endif
						enddo
					endif
				enddo
			enddo
		enddo
		
		beta = 0.d0
		alpha = 0.5d0 * alpha
						
		np = in_params%n_total
		call svdcmp(alpha,np,np,np,np,w,v)

! Invert the diagonal matrix
		wmax = maxval(w)
		do i = 1, np
			if (w(i) < wmax * 1.d-6) then
				beta(i,i) = 0.d0
			else
				beta(i,i) = 1.0 / w(i)
			endif
		enddo
		
		alpha = matmul(matmul(v, beta), transpose(alpha))
		
		do i = 1, np			
			if (params%inverted(i) == 1) then
				value = sqrt(abs(alpha(i,i)) * chi2Level)
				select case(i)
					case(1) 
						error%bgauss = value
					case(2) 
						error%thetabd = value
					case(3) 
						error%chibd = value
					case(4) 
						error%vdopp = value
					case(5) 
						error%dtau = value
					case(6) 
						error%delta_collision = value
					case(7)
						error%vmacro = value
					case(8)
						error%damping = value
					case(9)
						error%beta = value
					case(10)
						error%height = value
					case(11)
						error%dtau2 = value
					case(12)
						error%vmacro2 = value
					case(13) 
						error%bgauss2 = value
					case(14) 
						error%thetabd2 = value
					case(15) 
						error%chibd2 = value
					case(16) 
						error%vdopp2 = value
					case(17) 
						error%ff = value
					case(18)
						error%beta2 = value
				end select								
			endif			
		enddo
						
		deallocate(alpha)
		deallocate(beta)
		deallocate(v)
		deallocate(w)
		deallocate(x)
	
	end subroutine compute_uncertainty

	
!------------------------------------------------------------
! Re-scale the parameters to make them of the order of unity
!------------------------------------------------------------
	subroutine compress_parameters(in_params,in_params_scaled)
	type(variable_parameters) :: in_params, in_params_scaled
		in_params_scaled%bgauss = in_params%bgauss / 1000.d0
		in_params_scaled%thetabd = in_params%thetabd / 10.d0
		in_params_scaled%chibd = in_params%chibd / 10.d0
		in_params_scaled%vdopp = in_params%vdopp / 5.d0
		in_params_scaled%dtau = in_params%dtau / 1.5d0
		in_params_scaled%delta_collision = in_params%delta_collision
		in_params_scaled%vmacro = in_params%vmacro / 5.d0
		in_params_scaled%damping = in_params%damping
		in_params_scaled%beta = in_params%beta / 3.d0
		in_params_scaled%height = in_params%height / 30.d0

! Two components with same field and Doppler width
		if (in_params%nslabs /= 1) then
			in_params_scaled%dtau2 = in_params%dtau2 / 1.5d0
			in_params_scaled%vmacro2 = in_params%vmacro2 / 5.d0
		endif

! Two components with different field and Doppler width
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			in_params_scaled%bgauss2 = in_params%bgauss2 / 1000.d0
			in_params_scaled%thetabd2 = in_params%thetabd2 / 10.d0
			in_params_scaled%chibd2 = in_params%chibd2 / 10.d0
			in_params_scaled%vdopp2 = in_params%vdopp2 / 5.d0
		endif

! Two components in the same pixel
		if (in_params%nslabs == -2) then
			in_params_scaled%ff = in_params%ff
		endif

		if (in_params%nslabs == 3 .or. in_params%nslabs == -2 .or. in_params%nslabs == 2) then
			in_params_scaled%beta2 = in_params%beta2 / 3.d0
		endif

		
	end subroutine compress_parameters	
	
!------------------------------------------------------------
! Re-scale the parameters to their original values
!------------------------------------------------------------
	subroutine expand_parameters(in_params_scaled,in_params)
	type(variable_parameters) :: in_params, in_params_scaled
		in_params%bgauss = in_params_scaled%bgauss * 1000.d0
		in_params%thetabd = in_params_scaled%thetabd * 10.d0
		in_params%chibd = in_params_scaled%chibd * 10.d0
		in_params%vdopp = in_params_scaled%vdopp * 5.d0
		in_params%dtau = in_params_scaled%dtau * 1.5d0
		in_params%delta_collision = in_params_scaled%delta_collision
		in_params%vmacro = in_params_scaled%vmacro * 5.d0
		in_params%damping = in_params_scaled%damping
		in_params%beta = in_params_scaled%beta * 3.d0
		in_params%height = in_params_scaled%height * 30.d0

! Two components with same field and Doppler width
		if (in_params%nslabs /= 1) then
			in_params%dtau2 = in_params_scaled%dtau2 * 1.5d0
			in_params%vmacro2 = in_params_scaled%vmacro2 * 5.d0
		endif

! Two components with different field and Doppler width
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			in_params%bgauss2 = in_params_scaled%bgauss2 * 1000.d0
			in_params%thetabd2 = in_params_scaled%thetabd2 * 10.d0
			in_params%chibd2 = in_params_scaled%chibd2 * 10.d0
			in_params%vdopp2 = in_params_scaled%vdopp2 * 5.d0
		endif

! Two components in the same pixel
		if (in_params%nslabs == -2) then
			in_params%ff = in_params_scaled%ff
		endif

		if (in_params%nslabs == 3 .or. in_params%nslabs == -2 .or. in_params%nslabs == 2) then
			in_params%beta2 = in_params_scaled%beta2 * 3.d0
		endif

		
	end subroutine expand_parameters	



!*************************************************************
!*************************************************************
! DIRECT METHOD
!*************************************************************
!*************************************************************

!------------------------------------------------------------
! Invert some parameters with the DIRECT method
!------------------------------------------------------------
	subroutine invert_with_direct(out_params, in_fixed, myrank, synth_option, error)
	real(kind=8) :: DIReps, DIRf
	integer :: ndim, DIRmaxf, DIRmaxT, myrank, DIRmaxf_input, synth_option
	integer :: DIRalg, error
	integer :: IError, logfile
	real(kind=8) :: fglobal, fglper, volper, sigmaper, volper_input
	real(kind=8), allocatable :: u(:), l(:)
	integer :: n, i, j
	real(kind=8), allocatable :: DIRx(:)
	
	integer, parameter :: iisize=300, idsize=300, icsize=30
	integer :: iidata(iisize), verbose
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	character(len=120) :: location_file_direct
	
	type(variable_parameters) :: out_params
	type(fixed_parameters) :: in_fixed
						
! 		open(unit=23,file='logfile',action='write',status='replace')

		error = 0
		logfile = 23
		fglobal = -1.d100
		fglper = 0.d0
		sigmaper  = -1.d0
		
		DIReps = 1.d-6
		DIRmaxT = 1000
		
		ndim = sum(params%inverted)
		allocate(u(ndim))
		allocate(l(ndim))
		allocate(DIRx(ndim))
		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				l(j) = in_fixed%lower_direct(i)
				u(j) = in_fixed%upper_direct(i)
 				! print *, 'Inverting ', trim(adjustl(parameters_name(i))), l(j), u(j)
				j = j + 1
			endif
		enddo
				
		DIRalg = 1

		if (myrank == 0) then
			verbose = 1
		else
			verbose = 0
		endif

		iidata(1) = myrank

		DIRmaxf_input = in_fixed%DIRmaxf
		volper_input = in_fixed%volper

		if (synth_option == 2) then			
			call DIRECT(fcn,DIRx,ndim,DIReps,DIRmaxf_input,DIRmaxT, DIRf, l, u, DIRalg, Ierror, logfile, &
				fglobal, fglper, volper_input, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize, verbose)			
		else
			call DIRECT(fcn_simplified_StokesI,DIRx,ndim,DIReps,DIRmaxf_input,DIRmaxT, DIRf, l, u, DIRalg, Ierror, logfile, &
				fglobal, fglper, volper_input, sigmaper, iidata, iisize, ddata, idsize, cdata, icsize, verbose)
		endif
		
! If Ierror < 0, then some fatal error occured
		if (Ierror < 0) then
			error = 1
			return
		endif

! If Ierror > 0, then
		if (Ierror > 0) then
			if (Ierror == 1) then
!				write(*,*) 'DIRECT - Number of function evaluations done is larger than ', in_fixed%DIRmaxf
			endif
			if (Ierror == 2) then
!				write(*,*) 'DIRECT - Number of iterations is equal to ', DIRmaxT
			endif
			if (Ierror == 3) then
!				write(*,*) 'DIRECT - The best function value found is within ', fglper
			endif
			if (Ierror == 4) then
!				write(*,*) 'DIRECT - The volume of the hyperrectangle is less than ', in_fixed%volper
			endif
			if (Ierror == 5) then
!				write(*,*) 'DIRECT - The measure of the hyperrectangle is less than ', sigmaper
			endif
		endif
				
			
		j = 1		
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1) 
						out_params%bgauss = DIRx(j)
					case(2) 
						out_params%thetabd = DIRx(j)
					case(3) 
						out_params%chibd = DIRx(j)
					case(4) 
						out_params%vdopp = DIRx(j)
					case(5) 
						out_params%dtau = DIRx(j)
					case(6) 
						out_params%delta_collision = DIRx(j)
					case(7)
						out_params%vmacro = DIRx(j)
					case(8)
						out_params%damping = DIRx(j)
					case(9)
						out_params%beta = DIRx(j)
					case(10)
						out_params%height = DIRx(j)
					case(11)
						out_params%dtau2 = DIRx(j)
					case(12)
						out_params%vmacro2 = DIRx(j)
					case(13) 
						out_params%bgauss2 = DIRx(j)
					case(14) 
						out_params%thetabd2 = DIRx(j)
					case(15) 
						out_params%chibd2 = DIRx(j)
					case(16) 
						out_params%vdopp2 = DIRx(j)
					case(17)
						out_params%ff = DIRx(j)
					case(18) 
						out_params%beta2 = DIRx(j)
				end select
				j = j + 1
			endif
		enddo
			
! 		close(23)

		deallocate(u)
		deallocate(l)
		deallocate(DIRx)
			
	end subroutine invert_with_direct
	
!------------------------------------------------------------
! Function that returns the the DIRECT
!------------------------------------------------------------
	subroutine fcn(n, x, f, flag, iidata, iisize, ddata, idsize, cdata, icsize)
	integer :: n,flag
   real(kind=8) :: x(n)
	real(kind=8) :: f
	integer :: iisize, idsize, icsize, i, j, error
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	type(variable_parameters) :: trial
	character(len=120) :: temporal_file

		temporal_file = 'temporal.prof'

		
		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1) 
						params%bgauss = x(j)
					case(2) 
						params%thetabd = x(j)
					case(3) 
						params%chibd = x(j)
					case(4) 
						params%vdopp = x(j)
					case(5) 
						params%dtau = x(j)
					case(6) 
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
					case(13) 
						params%bgauss2 = x(j)
					case(14) 
						params%thetabd2 = x(j)
					case(15) 
						params%chibd2 = x(j)
					case(16) 
						params%vdopp2 = x(j)
					case(17)
						params%ff = x(j)
					case(18) 
						params%beta2 = x(j)
				end select
				j = j + 1
			endif
		enddo
		
			
		call do_synthesis(params, fixed, observation, inversion%stokes_unperturbed, error)

! If the computation is carried out by the master
		if (iidata(1) == 0) then
			call write_final_profiles(temporal_file, observation, inversion)
		endif
		
		f = compute_chisq(observation,inversion)
		
!   		write(*,FMT='(10X,A,I4,A,F18.8)') 'D  - chi^2(', iidata(1), ') : ', f
		
		if (iidata(1) == 0) then
 			call print_parameters(params,'      -Parameters : ',.TRUE.)
		
			print *, 'chi^2 : ', f
			write(25,FMT='(11(E15.5,2X))') x, f
		endif

		flag = 0

! If there was an error, return a non-zero flag
		if (error /= 0) flag = 1

	end subroutine fcn

!------------------------------------------------------------
! Function that returns the chi^2 for the DIRECT problem. This solves an approximate RT problem for Stokes I
!------------------------------------------------------------
	subroutine fcn_simplified_StokesI(n, x, f, flag, iidata, iisize, ddata, idsize, cdata, icsize)
	integer :: n,flag
   real(kind=8) :: x(n)
	real(kind=8) :: f
	integer :: iisize, idsize, icsize, i, j, error
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	type(variable_parameters) :: trial
	character(len=120) :: temporal_file
	real(kind=8), allocatable :: onum(:), prof(:), prof2(:)
	real(kind=8) :: va, dnum, adamp, onum0, onum1, onum2

		temporal_file = 'temporal.prof'

		j = 1
		do i = 1, params%n_total
			if (params%inverted(i) == 1) then
				select case(i)
					case(1)
						params%bgauss = x(j)
					case(2)
						params%thetabd = x(j)
					case(3)
						params%chibd = x(j)
					case(4)
						params%vdopp = x(j)
					case(5)
						params%dtau = x(j)
					case(6)
						params%delta_collision = x(j)
					case(7)
						params%vmacro = x(j)
					case(8)
						params%damping = x(j)
					case(9)
						params%beta = x(j)
					case(10)
						params%height = x(j)
					case(11)
						params%dtau2 = x(j)
					case(12)
						params%vmacro2 = x(j)
					case(13)
						params%bgauss2 = x(j)
					case(14)
						params%thetabd2 = x(j)
					case(15)
						params%chibd2 = x(j)
					case(16)
						params%vdopp2 = x(j)
					case(17)
						params%ff = x(j)
					case(18)
						params%beta2 = x(j)
				end select
				j = j + 1
			endif
		enddo

		allocate(onum(fixed%no))
		allocate(prof(fixed%no))

		onum = -1.d8 * observation%wl / fixed%wl**2

		inversion%stokes_unperturbed = 0.d0

! Only one slab
		if (params%nslabs == 1) then
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			inversion%stokes_unperturbed(0,:) = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = 1.d0 - exp(-params%dtau / 1.56d0 * inversion%stokes_unperturbed(0,:))
			else
				inversion%stokes_unperturbed(0,:) = exp(-params%dtau / 1.56d0 * inversion%stokes_unperturbed(0,:))
			endif			
		endif

! Two slabs
 		if (params%nslabs == 2 .or. params%nslabs == 3) then
! First slab
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = 1.d0 - exp(-params%dtau / 1.56d0 * prof)
			else
				inversion%stokes_unperturbed(0,:) = exp(-params%dtau / 1.56d0 * prof)
			endif

! Second slab
! If nslabs=3, then we use a different Doppler width. If not, we use the same for both components
			if (params%nslabs == 3) then
				dnum = params%vdopp2*1.d5 / (fixed%wl*1.d-8*PC)
				if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp2*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			endif
			va = params%vmacro2*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			inversion%stokes_unperturbed(0,:) = exp(-params%dtau2 / 1.56d0 * prof) * inversion%stokes_unperturbed(0,:)

 		endif

! Two slabs on the same pixel
		if (params%nslabs == -2) then

			allocate(prof2(fixed%no))

! First slab
			dnum = params%vdopp*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

! Second slab
			dnum = params%vdopp2*1.d5 / (fixed%wl*1.d-8*PC)
			va = params%vmacro2*1.d5 / (fixed%wl*1.d-8*PC)
			adamp = params%damping

 			if (fixed%damping_treatment == 1) then
 				adamp = fixed%wl * 1.d-8 / (params%vdopp2*1.d5) * aesto(fixed%nemiss) * abs(params%damping)
 			endif

			onum0 = -1.26d0*1.d-8 / (fixed%wl*1.d-8)**2
			onum1 = (-1.26d0+0.09d0)*1.d-8 / (fixed%wl*1.d-8)**2
			onum2 = 0.d0
			prof2 = profile(adamp,(onum0-onum-va)/dnum) + &
				0.6*profile(adamp,(onum1-onum-va)/dnum) + &
				0.2*profile(adamp,(onum2-onum-va)/dnum)

			if (fixed%Stokes_incident(0) == 0) then
				inversion%stokes_unperturbed(0,:) = params%ff * (1.d0 - exp(-params%dtau / 1.56d0 * prof)) + &
					(1.d0-params%ff) * (1.d0 - exp(-params%dtau2 / 1.56d0 * prof2))
			else
				inversion%stokes_unperturbed(0,:) = params%ff * exp(-params%dtau / 1.56d0 * prof) + &
					(1.d0-params%ff) * exp(-params%dtau2 / 1.56d0 * prof2)
			endif

			deallocate(prof2)
		endif

		f = compute_chisq(observation,inversion)

		deallocate(onum)
		deallocate(prof)

		if (iidata(1) == 0) then
			call write_final_profiles(temporal_file, observation, inversion)
		endif

! 		write(*,FMT='(10X,A,I4,A,F18.8)') 'D  - chi^2(', iidata(1), ') : ', f

		if (iidata(1) == 0) then
 			call print_parameters(params,'      -Parameters : ',.TRUE.)

			print *, 'chi^2 : ', f
			write(25,FMT='(11(E15.5,2X))') x, f
		endif

		flag = 0

! If there was an error, return a non-zero flag
! 		if (error /= 0) flag = 1

	end subroutine fcn_simplified_StokesI

end module marquardt
