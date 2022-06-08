module rt_coef
use maths
use vars
implicit none
contains

!------------------------------------------------------------
! Calculate the radiative transfer coefficients
!------------------------------------------------------------
	subroutine calc_rt_coef(in_params,in_fixed,in_observation,component)
	type(variable_parameters) :: in_params
	type(fixed_parameters) :: in_fixed
	type(type_observation) :: in_observation
	integer :: k, q, p, j, component
	real(kind=8) :: theta, chi, gamma, sign
	complex(kind=8) :: tr(0:2,-2:2,0:3), trp(0:2,-2:2,0:3), ii, suma
	real(kind=8) :: onum0, dnum, adamp, reduc(0:2,-2:2,-2:2)
	integer, allocatable :: njlevu(:), njlevl(:)
	real(kind=8), allocatable :: autl(:,:), autu(:,:), cl(:,:,:), cu(:,:,:), e0(:)
	real(kind=8), allocatable :: tmp1(:), tmp2(:), onum(:)
	integer :: nl, nu, ll2, lu2, jminl2, jmaxl2, jminu2, jmaxu2, j2, j2max, njlargest
	integer :: kmax, qmax, i, jp2, kmin, qp
	complex(kind=8), allocatable :: rot_mat_vert_mag(:,:,:), prof(:)
	real(kind=8) :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, rimax, pq, pu, pv, thb, chb, freq, bul, va, blu, prod, bfield
	integer :: nloop, ml2, mu2, q2, mup2, qp2, bigq2, bigq, bigqu2, bigqu, jal2, jl2, jlp2, jau2, ju2, jus2, jaup2, jup2
	integer :: jsmalll, jsmallu, k2, jumin, kumax, ku, ku2, istok, kumin, mlp2, bigql2, bigql, jls2, jalp2, klmin
	integer :: klmax, kl, kl2, nfreq, tnumber
	
	integer :: values_start(8), values_end(8)
	real(kind=8) :: start_time, end_time
	
	
		ii = cmplx(0.d0,1.d0)
		
      theta = in_fixed%thetad * PI/180.d0
      chi = in_fixed%chid * PI/180.d0
      gamma = in_fixed%gammad * PI/180.d0
		
! The emission and absorption formulae are given in the magnetic reference frame. Therefore, we need the rho^K_Q
! calculated in this reference frame. The T^K_Q(i,Omega) tensor thus gives us the emission and absorption coefficient
! when observing at a given direction Omega with respect with the magnetic field reference frame. Then, it is composed
! of two rotations: one to carry the observer reference frame to the vertical, given by the rotation R1=(-gamma,-theta,-chi)
! and the other one to carry this vertical reference frame to the magnetic reference frame, given by the 
! rotation (chb,thb,0)

! Which component we are computing
		if (component == 1) then
			thb = in_params%thetabd * PI / 180.d0
			chb = in_params%chibd * PI / 180.d0
			bfield = in_params%bgauss
		else
			thb = in_params%thetabd2 * PI / 180.d0
			chb = in_params%chibd2 * PI / 180.d0
			bfield = in_params%bgauss2
		endif
		
! Random azimuth solution. Since the density matrix in the magnetic field
! reference frame is independent on azimuth, in the random azimuth solution
! I think one can do the solution for any azimuth. The reason is that, even 
! if the SEE are solved for the density matrix in the vertical reference frame,
! they are later transformed to the magnetic field reference frame before 
! calculating the emission/absorption coefficients
		if (in_params%chibd == 999.d0) chb = 40.d0 * PI / 180.d0
		
      call tkq2(-gamma,-theta,-chi,chb,thb,0.d0,tr)
      
! Random azimuth solution. It follows the descriptions of 
! Belluzzi, Trujillo Bueno & Landi Degl'Innocenti ApJ 666, 588 (2007)
      if (in_params%chibd == 999.d0) then
      	print *, 'Random azimuth solution...'
      	call reduced_matrix(2,thb,reduc)
      	call tkq2(-gamma,-theta,-chi,0.d0,0.d0,0.d0,trp)
      	do k = 0, 2
      		do bigq = -2, 2
      			do istok = 0, 3
      				tr(k,bigq,istok) = trp(k,0,istok) * reduc(k,0,bigq)
      			enddo
      		enddo
      	enddo
      endif

		if (.not.associated(in_fixed%epsilon)) allocate(in_fixed%epsilon(0:3,in_fixed%no))
		if (.not.associated(in_fixed%epsilon_zeeman)) allocate(in_fixed%epsilon_zeeman(0:3,in_fixed%no))
		if (.not.associated(in_fixed%eta)) allocate(in_fixed%eta(0:3,in_fixed%no))
		if (.not.associated(in_fixed%eta_zeeman)) allocate(in_fixed%eta_zeeman(0:3,in_fixed%no))
		if (.not.associated(in_fixed%eta_stim)) allocate(in_fixed%eta_stim(0:3,in_fixed%no))
		if (.not.associated(in_fixed%eta_stim_zeeman)) allocate(in_fixed%eta_stim_zeeman(0:3,in_fixed%no))
		if (.not.associated(in_fixed%mag_opt)) allocate(in_fixed%mag_opt(0:3,in_fixed%no))
		if (.not.associated(in_fixed%mag_opt_zeeman)) allocate(in_fixed%mag_opt_zeeman(0:3,in_fixed%no))
		if (.not.associated(in_fixed%mag_opt_stim)) allocate(in_fixed%mag_opt_stim(0:3,in_fixed%no))
		if (.not.associated(in_fixed%mag_opt_stim_zeeman)) allocate(in_fixed%mag_opt_stim_zeeman(0:3,in_fixed%no))
		
		in_fixed%epsilon = 0.d0
		in_fixed%epsilon_zeeman = 0.d0
		in_fixed%eta = 0.d0
		in_fixed%eta_zeeman = 0.d0
		in_fixed%eta_stim = 0.d0
		in_fixed%eta_stim_zeeman = 0.d0
		
		in_fixed%mag_opt = 0.d0
		in_fixed%mag_opt_stim = 0.d0
		in_fixed%mag_opt_zeeman = 0.d0
		in_fixed%mag_opt_stim_zeeman = 0.d0
		if (.not.allocated(onum)) allocate(onum(in_fixed%no))
		
!EDGAR: wl should be here wavelength, but these expressions reminds me those in wavenumber: REMEMBER CHECK
! In the synthesis mode, generate a new wavelength axis (wavenumber in this case)
! #if ! defined(python)
! 		if (working_mode == 0) then		
! 			do i = 1, in_fixed%no
! 				onum(i) = in_fixed%omin + (in_fixed%omax-in_fixed%omin) * dfloat(i-1) / dfloat(in_fixed%no-1)
! 			enddo			
! 			in_observation%wl = -1.d-8 * in_fixed%wl**2 * onum
! 		endif
! #else
		onum = -1.d8 * in_observation%wl / in_fixed%wl**2
! #endif

! In the inversion mode, use the wavelength axis (wavenumber in this case) from the observation
		if (working_mode == 1) then
			onum = -1.d8 * in_observation%wl / in_fixed%wl**2
		endif
				
! Transform the Doppler velocity in wavenumber and we calculate the reduced damping constant

! Which component we are computing
		if (component == 1) then
			dnum = in_params%vdopp*1.d5 / (in_fixed%wl*1.d-8*PC)		! Delta w = v_th / (lambda*c) (from eqs. 5.43 of the book)
		else
			dnum = in_params%vdopp2*1.d5 / (in_fixed%wl*1.d-8*PC)		! Delta w = v_th / (lambda*c) (from eqs. 5.43 of the book)
		endif
		
		adamp = in_params%damping

! If we compute the damping parameter using the natural broadening, use Eq. 11.47 and 11.15 of Grey "Stellar photospheres"
! Modify the number with an enhancement given by the absolute value of the (negative) damping
! The natural damping is given by gamma=4*pi*Aul and the damping constant is a=gamma / (4*pi) * lambda0^2/c * 1/dlambdad
		if (in_fixed%damping_treatment == 1) then
			adamp = in_fixed%wl * 1.d-8 / (in_params%vdopp*1.d5) * aesto(in_fixed%nemiss) * abs(in_params%damping)
		endif

		if (component == 1) then
			va = in_params%vmacro * 1.d5 / (in_fixed%wl*1.d-8*PC)
		else
			va = in_params%vmacro2 * 1.d5 / (in_fixed%wl*1.d-8*PC)
		endif
				
!-------------------------------------------------------------------------
!----- We calculate eigenvalues and eigenvectors
!-------------------------------------------------------------------------
      nl = ntlsto(in_fixed%nemiss)
      nu = ntusto(in_fixed%nemiss)
      ll2 = lsto2(nl)
      lu2 = lsto2(nu)
      jminl2 = abs(ll2-is2)
      jmaxl2 = ll2+is2
      jminu2 = iabs(lu2-is2)
      jmaxu2 = lu2+is2
				
! Introduce the energy of the levels of the lower level
		allocate(e0(0:jlimit2))
      do j2 = jminl2, jmaxl2, 2
      	e0(j2) = energy(nl,j2)			
		enddo
		
! Calculate the size of the eigenvalue and eigenvectors array
		call paschen_size(ll2,is2,j2max,njlargest)
						
		if (.not.allocated(njlevl)) allocate(njlevl(-j2max:j2max))
		if (.not.allocated(autl)) allocate(autl(-j2max:j2max,njlargest))
		autl = 0.d0
		if (.not.allocated(cl)) allocate(cl(-j2max:j2max,njlargest,0:j2max))
		cl = 0.d0		

! Diagonalize the Paschen-Back Hamiltonian
      call paschen(ll2,is2,e0,bfield,njlevl,autl,cl,j2max,njlargest)
            		
! And now the same for the upper level
      do j2 = jminu2, jmaxu2, 2
      	e0(j2) = energy(nu,j2)
		enddo
		
		call paschen_size(lu2,is2,j2max,njlargest)
		
		if (.not.allocated(njlevu)) allocate(njlevu(-j2max:j2max))
		if (.not.allocated(autu)) allocate(autu(-j2max:j2max,njlargest))
		autu = 0.d0
		if (.not.allocated(cu)) allocate(cu(-j2max:j2max,njlargest,0:j2max))
		cu = 0.d0
		
      call paschen(lu2,is2,e0,bfield,njlevu,autu,cu,j2max,njlargest)
		
		deallocate(e0)
		
! Write a matrix with the value of the rho^K_Q(J,J') of the upper and lower term
		j2max = maxval(j2tab)
		kmax = maxval(ktab)
		qmax = maxval(qtab)
		if (.not.allocated(rhou)) allocate(rhou(0:j2max,0:j2max,0:kmax,-qmax:qmax))
		rhou = 0.d0
		if (.not.allocated(rhomu)) allocate(rhomu(0:j2max,0:j2max,0:kmax,-qmax:qmax))
		rhomu = 0.d0
		if (.not.allocated(rhol)) allocate(rhol(0:j2max,0:j2max,0:kmax,-qmax:qmax))
		rhol = 0.d0
		if (.not.allocated(rhoml)) allocate(rhoml(0:j2max,0:j2max,0:kmax,-qmax:qmax))
		rhoml = 0.d0
		
		do i = 1, nrhos
			if (ntab(i) == nu) then
				if (irtab(i) == 1) then
					call set_real(rhou(j2tab(i),jp2tab(i),ktab(i),qtab(i)),SEE_b(i))
				else
					call set_imaginary(rhou(j2tab(i),jp2tab(i),ktab(i),qtab(i)),SEE_b(i))
				endif
				
				if (j2tab(i) == jp2tab(i)) then
          		if (qtab(i) == 0) then
              		call set_imaginary(rhou(j2tab(i),jp2tab(i),ktab(i),qtab(i)),0.d0)
              	else
						sign = 1.d0
              		if (mod(qtab(i),2) /= 0) sign=-1.d0              		
						if(irtab(i) == 1) then
							call set_real(rhou(j2tab(i),jp2tab(i),ktab(i),-qtab(i)),sign*SEE_b(i))
						else
							call set_imaginary(rhou(j2tab(i),jp2tab(i),ktab(i),-qtab(i)),-sign*SEE_b(i))
						endif
					endif
				else
      			sign = 1.d0				
      			if (mod(j2tab(i)-jp2tab(i)-2*qtab(i),4) /= 0) sign = -1.d0
      			if (irtab(i) == 1) then				
						call set_real(rhou(jp2tab(i),j2tab(i),ktab(i),-qtab(i)),sign*SEE_b(i))
      			else
						call set_imaginary(rhou(jp2tab(i),j2tab(i),ktab(i),-qtab(i)),-sign*SEE_b(i))
					endif
				endif
			endif
			if (ntab(i) == nl) then
				if (irtab(i) == 1) then
					call set_real(rhol(j2tab(i),jp2tab(i),ktab(i),qtab(i)),SEE_b(i))
				else
					call set_imaginary(rhol(j2tab(i),jp2tab(i),ktab(i),qtab(i)),SEE_b(i))
				endif
				
				if (j2tab(i) == jp2tab(i)) then
          		if (qtab(i) == 0) then
              		call set_imaginary(rhol(j2tab(i),jp2tab(i),ktab(i),qtab(i)),0.d0)
              	else
						sign = 1.d0
              		if (mod(qtab(i),2) /= 0) sign=-1.d0              		
						if(irtab(i) == 1) then
							call set_real(rhol(j2tab(i),jp2tab(i),ktab(i),-qtab(i)),sign*SEE_b(i))
						else
							call set_imaginary(rhol(j2tab(i),jp2tab(i),ktab(i),-qtab(i)),-sign*SEE_b(i))
						endif
					endif
				else
      			sign = 1.d0				
      			if (mod(j2tab(i)-jp2tab(i)-2*qtab(i),4) /= 0) sign = -1.d0
      			if (irtab(i) == 1) then				
						call set_real(rhol(jp2tab(i),j2tab(i),ktab(i),-qtab(i)),sign*SEE_b(i))
      			else
						call set_imaginary(rhol(jp2tab(i),j2tab(i),ktab(i),-qtab(i)),-sign*SEE_b(i))
					endif
				endif
			endif
		enddo
		
! We write the upper level rho^K_Q(J,J') in the reference frame of the vertical
! 		open(unit=12, file=output_rho_vertical_upper,action='write',status='replace')
! 		do j2 = jminu2, jmaxu2, 2
! 			do jp2 = jminu2, jmaxu2, 2
! 				kmin = abs(jp2-j2)/2
! 				kmax = (jp2+j2) / 2
! 				do k = kmin, kmax
! 					do q = -k, k
! 						write(12,FMT='(4(1X,I3),2(1x,e15.7))') j2, jp2, k, q, real(rhou(j2,jp2,k,q)), aimag(rhou(j2,jp2,k,q))
! 					enddo
! 				enddo
! 			enddo
! 		enddo
! 		close(12)
		
! We write the lower level rho^K_Q(J,J') in the reference frame of the vertical
! 		open(unit=12, file=output_rho_vertical_lower,action='write',status='replace')
! 		do j2 = jminl2, jmaxl2, 2
! 			do jp2 = jminl2, jmaxl2, 2
! 				kmin = abs(jp2-j2)/2
! 				kmax = (jp2+j2) / 2
! 				do k = kmin, kmax
! 					do q = -k, k
! 						write(12,FMT='(4(1X,I3),2(1x,e15.7))') j2, jp2, k, q, real(rhol(j2,jp2,k,q)), aimag(rhol(j2,jp2,k,q))
! 					enddo
! 				enddo
! 			enddo
! 		enddo
! 		close(12)		
		
! We transform the density matrix to the reference system of the magnetic field. We perform a rotation of the 
! density matrix elements following eq. 3.98 of Landi & Landolfi (2004)
! We need the rotation matrix D^K_{Q'Q}(R)^* = D^K_{QQ'}(-R)
		kmax = max(jmaxu2,jmaxl2)
		
		if (.not.allocated(rot_mat_vert_mag)) allocate(rot_mat_vert_mag(0:kmax,-kmax:kmax,-kmax:kmax))
		rot_mat_vert_mag = 0.d0

		call rotmatsu(kmax,0.d0,-thb,-chb,rot_mat_vert_mag)
		
		do j2 = jminu2, jmaxu2, 2
			do jp2 = jminu2, jmaxu2, 2
				kmin = abs(j2-jp2)/2
				kmax = (j2+jp2) / 2
				do k = kmin, kmax
					do q = -k, k
						suma = 0.d0
						do qp = -k, k
							suma = suma + rhou(j2,jp2,k,qp) * rot_mat_vert_mag(k,q,qp)							
						enddo
						rhomu(j2,jp2,k,q) = suma
					enddo
				enddo
			enddo
		enddo		
				
		do j2 = jminl2, jmaxl2, 2
			do jp2 = jminl2, jmaxl2, 2
				kmin = abs(j2-jp2)/2
				kmax = (j2+jp2) / 2
				do k = kmin, kmax
					do q = -k, k
						suma = 0.d0
						do qp = -k, k
							suma = suma + rhol(j2,jp2,k,qp) * rot_mat_vert_mag(k,q,qp)							
						enddo
						rhoml(j2,jp2,k,q) = suma						
					enddo
				enddo
			enddo
		enddo
		
! We write the upper level rho^K_Q(J,J') in the reference frame of the vertical
! 		open(unit=12, file=output_rho_magnetic_upper,action='write',status='replace')
! 		do j2 = jminu2, jmaxu2, 2
! 			do jp2 = jminu2, jmaxu2, 2
! 				kmin = abs(jp2-j2)/2
! 				kmax = (jp2+j2) / 2
! 				do k = kmin, kmax
! 					do q = -k, k
! 						write(12,FMT='(4(1X,I3),2(1x,e15.7))') j2, jp2, k, q, real(rhomu(j2,jp2,k,q)), aimag(rhomu(j2,jp2,k,q))
! 					enddo
! 				enddo
! 			enddo
! 		enddo
! 		close(12)
		
! We write the lower level rho^K_Q(J,J') in the reference frame of the vertical
! 		open(unit=12, file=output_rho_magnetic_lower,action='write',status='replace')
! 		do j2 = jminl2, jmaxl2, 2
! 			do jp2 = jminl2, jmaxl2, 2
! 				kmin = abs(jp2-j2)/2
! 				kmax = (jp2+j2) / 2
! 				do k = kmin, kmax
! 					do q = -k, k
! 						write(12,FMT='(4(1X,I3),2(1x,e15.7))') j2, jp2, k, q, real(rhoml(j2,jp2,k,q)), aimag(rhoml(j2,jp2,k,q))
! 					enddo
! 				enddo
! 			enddo
! 		enddo
! 		close(12)		

!-------------------------------------------------------------------------
! And we perform the summation given by eq. 7.47 of Landi & Landolfi (2004)
! Here we calculate only 7.47e (equivalent to 7.47b of stimulated emission but multiplied by 2*h*nu^3/c^2) 
! because we are interested in the emission coefficient
!-------------------------------------------------------------------------
!		call date_and_time(values=values_start)		
		nloop = 0
      x0 = dfloat(lu2+1)
		nfreq = in_fixed%no
		
		if (.not.allocated(tmp1)) allocate(tmp1(nfreq))
		if (.not.allocated(tmp2)) allocate(tmp2(nfreq))
		if (.not.allocated(prof)) allocate(prof(nfreq))

      do ml2 = -jmaxl2, jmaxl2, 2
      	do mu2 = -jmaxu2, jmaxu2, 2
	      	q2 = ml2-mu2
	      	if (iabs(q2) <= 2) then
		      	do mup2 = -jmaxu2, jmaxu2, 2
			      	qp2 = ml2-mup2
			      	if (iabs(qp2) <= 2) then
				      	bigq2 = q2-qp2
				      	bigq = bigq2/2
				      	bigqu2 = mup2-mu2
				      	bigqu = bigqu2/2
				      	jal2 = iabs(ml2)
				      	if (jal2 < jminl2) jal2 = jminl2
				      	do jl2 = jal2, jmaxl2, 2
					      	do jlp2 = jal2, jmaxl2, 2
						      	jau2=iabs(mu2)
						      	if (jau2 < jminu2) jau2 = jminu2
						      	do ju2 = jau2, jmaxu2, 2
							      	if (iabs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
										else
								      	do jus2 = jau2, jmaxu2, 2
									      	jaup2=iabs(mup2)
									      	if (jaup2 < jminu2) jaup2 = jminu2
									      	do jup2 = jaup2, jmaxu2, 2
									      		if (iabs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
													else
										      		x1 = 1.d0
										      		if (mod(2+jup2-mu2+qp2,4) /= 0) x1 = -1.d0
										      		x2 = dsqrt(dfloat((jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)))
										      		x3 = w3js(ju2,jl2,2,-mu2,ml2,-q2)
										      		x4 = w3js(jup2,jlp2,2,-mup2,ml2,-qp2)
										      		x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
										      		x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
										      		do jsmalll = 1, njlevl(ml2)
											      		do jsmallu = 1, njlevu(mu2)
												      		onum0 = autu(mu2,jsmallu) - autl(ml2,jsmalll)
! Evaluate the profile
																prof = profile(adamp,(onum0-onum-va)/dnum)																

												      		x7=cl(ml2,jsmalll,jl2)*cl(ml2,jsmalll,jlp2)*&
																	cu(mu2,jsmallu,ju2)*cu(mu2,jsmallu,jus2)
												      		do k = 0, 2
													      		k2 = k*2
													      		x8 = w3js(2,2,k2,q2,-qp2,-bigq2)
													      		kumin = iabs(jup2-jus2)/2
													      		kumax = (jup2+jus2)/2
													      		do ku = kumin, kumax
														      		ku2 = ku*2
														      		x9 = w3js(jup2,jus2,ku2,mup2,-mu2,-bigqu2)
														      		x10 = dsqrt(dfloat(3*(k2+1)*(ku2+1)))
																		
																		prod = x0*x1*x2*x3*x4*x5*x6*x7*x8*x9*x10
																		
! Take the real part because we are calculating emission/absorption coefficients
																		do istok = 0, 3																			
																			tmp1 = prod * &
																				real(tr(k,bigq,istok)*rhomu(jup2,jus2,ku,bigqu)*prof)																			
																			in_fixed%epsilon(istok,:) = in_fixed%epsilon(istok,:) + tmp1
																			
! Now the magneto-optical coefficients																			
																			tmp2 = prod * &
																				aimag(tr(k,bigq,istok)*rhomu(jup2,jus2,ku,bigqu)*prof)
																			in_fixed%mag_opt_stim(istok,:) = in_fixed%mag_opt_stim(istok,:) + tmp2

! Only the Zeeman contribution, produced by the rho^0_0 multipole moments and the contribution of the splitting
																			if (ku == 0 .and. bigqu == 0) then
																				in_fixed%epsilon_zeeman(istok,:) = in_fixed%epsilon_zeeman(istok,:) + tmp1
																			
! Only the Zeeman contribution, produced by the rho^0_0 multipole moments and the contribution of the splitting																			
																				in_fixed%mag_opt_stim_zeeman(istok,:) = in_fixed%mag_opt_stim_zeeman(istok,:) + tmp2
																			endif																			
																			nloop = nloop + 1																			
																		enddo ! Stokes parameter																		
																	enddo  ! Ku
																enddo  ! K
															enddo  ! ju
														enddo  ! jl
													endif
												enddo  ! Ju'
											enddo  ! Ju''
										endif
									enddo  ! Ju
								enddo  ! Jl'
							enddo  ! Jl
						endif
					enddo  ! Mu'
				endif
			enddo  ! Mu
		enddo  ! Ml		
		
! 		deallocate(tmp1)
! 		deallocate(tmp2)
! 		deallocate(prof)
		
!		call date_and_time(values=values_end)
!		start_time = values_start(5) * 3600 + values_start(6) * 60 + values_start(7) + 0.001 * values_start(8)
!		end_time = values_end(5) * 3600 + values_end(6) * 60 + values_end(7) + 0.001 * values_end(8)
!		print *, 'Time1 : ', end_time-start_time, nloop
		
!		call date_and_time(values=values_start)
		
!-------------------------------------------------------------------------
! And we perform the summation given by eq. 7.47 of Landi & Landolfi (2004)
! Here we calculate only 7.47a because we are interested in the absorption coefficient
!-------------------------------------------------------------------------
		if (synthesis_mode /= 0) then
		
			nloop = 0
      	x0 = dfloat(lu2+1)   ! It is (2Lu+1) because we need B_lu, which is equal to (2Lu+1)/(2Ll+1)*B_ul

			if (.not.allocated(tmp1)) allocate(tmp1(nfreq))
			if (.not.allocated(tmp2)) allocate(tmp2(nfreq))
			if (.not.allocated(prof)) allocate(prof(nfreq))

      	do ml2 = -jmaxl2, jmaxl2, 2
      		do mu2 = -jmaxu2, jmaxu2, 2
	      		q2 = ml2-mu2
	      		if (iabs(q2) <= 2) then
		      		do mlp2 = -jmaxl2, jmaxl2, 2
			      		qp2 = mlp2-mu2
			      		if (iabs(qp2) <= 2) then
				      		bigq2 = q2-qp2
				      		bigq = bigq2/2
				      		bigql2 = ml2-mlp2
				      		bigql = bigql2/2
				      		jau2 = iabs(mu2)
				      		if (jau2 < jminu2) jau2 = jminu2
				      		do ju2 = jau2, jmaxu2, 2
					      		do jup2 = jau2, jmaxu2, 2
						      		jal2=iabs(ml2)
						      		if (jal2 < jminl2) jal2 = jminl2
						      		do jl2 = jal2, jmaxl2, 2
							      		if (iabs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
											else
								      		do jls2 = jal2, jmaxl2, 2
									      		jalp2=iabs(mlp2)
									      		if (jalp2 < jminl2) jalp2 = jminl2
									      		do jlp2 = jalp2, jmaxl2, 2
									      			if (iabs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
														else
										      			x1 = 1.d0
										      			if (mod(2+jls2-ml2+qp2,4) /= 0) x1 = -1.d0
										      			x2 = dsqrt(dfloat((jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)))
										      			x3 = w3js(ju2,jl2,2,-mu2,ml2,-q2)
										      			x4 = w3js(jup2,jlp2,2,-mu2,mlp2,-qp2)
										      			x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
										      			x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
										      			do jsmalll = 1, njlevl(ml2)
											      			do jsmallu = 1, njlevu(mu2)
												      			onum0 = autu(mu2,jsmallu) - autl(ml2,jsmalll)
	! Evaluate the profile
																	prof = profile(adamp,(onum0-onum-va)/dnum)

												      			x7=cl(ml2,jsmalll,jl2)*cl(ml2,jsmalll,jls2)*&
																		cu(mu2,jsmallu,ju2)*cu(mu2,jsmallu,jup2)
												      			do k = 0, 2
													      			k2 = k*2
													      			x8 = w3js(2,2,k2,q2,-qp2,-bigq2)
													      			klmin = iabs(jlp2-jls2)/2
													      			klmax = (jlp2+jls2)/2
													      			do kl = klmin, klmax
														      			kl2 = kl*2
														      			x9 = w3js(jls2,jlp2,kl2,ml2,-mlp2,-bigql2)
														      			x10 = dsqrt(dfloat(3*(k2+1)*(kl2+1)))

																			prod = x0*x1*x2*x3*x4*x5*x6*x7*x8*x9*x10
																			
	! Take the real part because we are calculating emission/absorption coefficients
																			do istok = 0, 3
																				tmp1 = prod * &
																					real(tr(k,bigq,istok)*rhoml(jls2,jlp2,kl,bigql)*prof)
																				in_fixed%eta(istok,:) = in_fixed%eta(istok,:) + tmp1

	! Now the magneto-optical coefficients
																				tmp2 = prod * &
																					aimag(tr(k,bigq,istok)*rhoml(jls2,jlp2,kl,bigql)*prof)
																				in_fixed%mag_opt(istok,:) = in_fixed%mag_opt(istok,:) + tmp2

	! Only the Zeeman contribution, produced by the rho^0_0 multipole moments and the contribution of the splitting
																				if (kl == 0 .and. bigql == 0) then
																					in_fixed%eta_zeeman(istok,:) = in_fixed%eta_zeeman(istok,:) + tmp1

	! Only the Zeeman contribution, produced by the rho^0_0 multipole moments and the contribution of the splitting
																					in_fixed%mag_opt_zeeman(istok,:) = in_fixed%mag_opt_zeeman(istok,:) + tmp2
																				endif
																				nloop = nloop + 1
																			enddo ! Stokes parameter											
																		enddo  ! Kl
																	enddo  ! K
																enddo  ! ju
															enddo  ! jl
														endif
													enddo  ! Jl'
												enddo  ! Jl''
											endif
										enddo  ! Jl
									enddo  ! Ju'
								enddo  ! Ju
							endif
						enddo  ! Ml'
					endif
				enddo  ! Mu
			enddo  ! Ml		
			
! 			deallocate(tmp1)
! 			deallocate(tmp2)
! 			deallocate(prof)						
		
		endif
		
!		call date_and_time(values=values_end)
!		start_time = values_start(5) * 3600 + values_start(6) * 60 + values_start(7) + 0.001 * values_start(8)
!		end_time = values_end(5) * 3600 + values_end(6) * 60 + values_end(7) + 0.001 * values_end(8)
!		print *, 'Time2 : ', end_time-start_time, nloop

! Introduce the hnu/4Pi factor
		freq = PC / (1.d-8*in_fixed%wl)
		Bul = PC**2 / (2.d0*PH*freq**3) * aesto(in_fixed%nemiss)
		Blu = (lu2+1.d0) / (ll2+1.d0) * Bul
		
! Save the value of the Einstein coefficients		
		in_fixed%Bul = Bul
		in_fixed%Blu = Blu
		in_fixed%Aul = aesto(in_fixed%nemiss)
		in_fixed%nu = freq
		
! Stimulated emission
!		eta_stim = PH * freq / (4.d0*PI) * Bul * epsilon
!		eta_stim_zeeman = PH * freq / (4.d0*PI) * Bul * epsilon_zeeman		
		
		in_fixed%eta_stim = in_fixed%epsilon
		in_fixed%eta_stim_zeeman = in_fixed%epsilon_zeeman		

! Emissivity
!		epsilon = PH * freq / (4.d0*PI) * Bul * (2.d0*PH*freq**3) / PC**2 * epsilon
!		epsilon_zeeman = PH * freq / (4.d0*PI) * Bul * (2.d0*PH*freq**3) / PC**2 * epsilon_zeeman
				
		in_fixed%epsilon = (2.d0*PH*freq**3) / PC**2 * in_fixed%epsilon
		in_fixed%epsilon_zeeman = (2.d0*PH*freq**3) / PC**2 * in_fixed%epsilon_zeeman

! Absorption
!		eta = PH * freq / (4.d0*PI) * Bul * eta
!		eta_zeeman = PH * freq / (4.d0*PI) * Bul * eta_zeeman		
				
		! eta = eta
		! eta_zeeman = eta_zeeman		
		
! Magneto-optical absorption
!		mag_opt = PH * freq / (4.d0*PI) * Bul * mag_opt
!		mag_opt_zeeman = PH * freq / (4.d0*PI) * Bul * mag_opt_zeeman
		
		! mag_opt = mag_opt
		! mag_opt_zeeman = mag_opt_zeeman
		
! Magneto-optical stimulated emission
!		mag_opt_stim = PH * freq / (4.d0*PI) * Bul * mag_opt_stim
!		mag_opt_stim_zeeman = PH * freq / (4.d0*PI) * Bul * mag_opt_stim_zeeman
		
		! mag_opt_stim = mag_opt_stim
		! mag_opt_stim_zeeman = mag_opt_stim_zeeman
		
!-------------------------------------------------------------------------		
! And finally write the emission coefficients		
!-------------------------------------------------------------------------				
! 		open(unit=12,file=output_rtcoef,action='write',status='replace')
! 		write(12,FMT='(A,56X,A,48X,A,48X,A,48X,A,48X,A,48X,A,48X)') 'Lambda','epsilon','eta_A','eta_S','rho_A','rho_S'
! 		do i = 1, in_fixed%no			
! 			write(12,FMT='(2(1X,E15.8),4(1X,E12.5),A2,4(1X,E12.5),A2,4(1X,E12.5),A2,3(1X,E12.5),A2,3(1X,E12.5))') &
! 				-1.d-8 * in_fixed%wl**2 * onum(i), &
! 				-1.d-8 * in_fixed%wl**2 * onum(i) + in_fixed%wl, (epsilon(j,i),j=0,3), '  ',(eta(j,i),j=0,3), '  ',&
! 				(eta_stim(j,i),j=0,3),'  ', (mag_opt(j,i),j=1,3), '  ', (mag_opt_stim(j,i),j=1,3)
! 		enddo
! 		close(12)
! 		
! 		open(unit=12,file=output_rtcoef_zeeman,action='write',status='replace')
! 		write(12,FMT='(A,56X,A,48X,A,48X,A,48X,A,48X,A,48X,A,48X)') 'Lambda','epsilon','eta_A','eta_S','rho_A','rho_S'
! 		do i = 1, in_fixed%no			
! 			write(12,FMT='(2(1X,E15.8),4(1X,E12.5),2X,4(1X,E12.5),2X,4(1X,E12.5),2X,3(1X,E12.5),2X,3(1X,E12.5))') &
! 				-1.d-8 * in_fixed%wl**2 * onum(i), &
! 				-1.d-8 * in_fixed%wl**2 * onum(i) + in_fixed%wl, (epsilon_zeeman(j,i),j=0,3), (eta_zeeman(j,i),j=0,3),&
! 				(eta_stim_zeeman(j,i),j=0,3), (mag_opt_zeeman(j,i),j=1,3), (mag_opt_stim_zeeman(j,i),j=1,3)
! 		enddo
! 		close(12)
		
		
! Clean the allocated memory
! #if defined(python)
		deallocate(tmp1)
		deallocate(tmp2)
		deallocate(onum)
		deallocate(prof)
		deallocate(njlevl)
		deallocate(autl)
		deallocate(cl)
		deallocate(njlevu)
		deallocate(autu)
		deallocate(cu)
		deallocate(rhou)
		deallocate(rhomu)
		deallocate(rhol)
		deallocate(rhoml)
		deallocate(rot_mat_vert_mag)
! #endif				
								
	end subroutine calc_rt_coef
end module rt_coef
