module io
! use vars
use maths
implicit none
contains

!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
	subroutine read_config
	character(len=100) :: config_file
	integer :: nargs, ios
		
		input_model_file = 'ATOMS/helium.mod'
		input_experiment = 'init_parameters.dat'		
		verbose_mode = 0		
		linear_solver = 0		
		synthesis_mode = 5		
		working_mode = 1
		
		
	end subroutine read_config
	

!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
	subroutine read_experiment(file,in_params,in_fixed)
	character(len=120) :: file
	type(variable_parameters) :: in_params
	type(fixed_parameters) :: in_fixed
	real(kind=8) :: delta, height
	integer :: j
	
		open(unit=12,file=file,action='read',status='old')

! Read if we want to include the stimulated emission
		call lb(12,5)
		read(12,*) isti
		if (isti == 1 .and. verbose_mode == 1) then
			print *, 'Including stimulated emission...'
		endif
	
! Read if we want to include the magnetic field
		call lb(12,2)
		read(12,*) imag
		if (imag == 1 .and. verbose_mode == 1) then
			print *, 'Including magnetic field ...'
		endif
		
! Read if we want to include the depolarizing collisions
		call lb(12,2)
		read(12,*) idep
		if (idep == 1 .and. verbose_mode == 1) then
			print *, 'Including depolarizing collisions ...'
		endif

! Value of the delta parameter for the depolarizing collisions		
		call lb(12,2)
		read(12,*) in_params%delta_collision
		
! Use Paschen-Back calculation or linear Zeeman effect
		call lb(12,2)
		read(12,*) use_paschen_back
		if (verbose_mode == 1) then
			if (use_paschen_back == 1) then
				print *, 'Taking into account the Paschen-Back effect...'
			else
				print *, 'Not taking into account the Paschen-Back effect...'
			endif
		endif

! Number of slabs
		call lb(12,2)
		read(12,*) in_params%nslabs

! Values characterizing the magnetic field
		call lb(12,2)
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			read(12,*) in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%bgauss2, in_params%thetabd2, in_params%chibd2
		else
			read(12,*) in_params%bgauss, in_params%thetabd, in_params%chibd
		endif
		
		if (imag == 1 .and. verbose_mode == 1) then
			print *, 'Field properties'
			write(*,FMT='(4X,A,F8.3,3X,A,F6.2,3X,A,F6.2)') 'B : ', in_params%bgauss, 'thetab : ', in_params%thetabd, 'chib : ',&
				in_params%chibd
		endif
		
! Values of the height for the calculation of the radiation anisotropy
		call lb(12,2)
		read(12,*) in_params%height
				
		if (verbose_mode == 1) then
			write(*,FMT='(A,I1,A)') 'Using ', abs(in_params%nslabs), ' slabs'
			if (in_params%nslabs == -2) then
				write(*,FMT='(A)') 'with filling factor'
			endif
		endif
		
! Value of the optical depth in the maximum (SLAB) or the strength of the line (ME)
		call lb(12,2)
		if (in_params%nslabs == 1) then
			read(12,*) in_params%dtau
		endif

		if (in_params%nslabs == 2 .or. in_params%nslabs == 3) then
			read(12,*) in_params%dtau, in_params%dtau2
		endif

! Two components with filling factor
		if (in_params%nslabs == -2) then
			read(12,*) in_params%dtau, in_params%dtau2, in_params%ff
		endif
		
! Value of the gradient of the source function (ME)
		call lb(12,2)
		read(12,*) in_params%beta
		
! Incident Stokes parameters
		call lb(12,2)
		read(12,*) (in_fixed%Stokes_incident(j),j=0,3)
		
! Transition where to compute the emission
		call lb(12,2)
		read(12,*) in_fixed%nemiss
		if (verbose_mode == 1) then			
			write(*,FMT='(A,3X,I3)') ' Calculating emission from line : ', in_fixed%nemiss
		endif
		
! Use atomic polarization
		call lb(12,2)
		read(12,*) in_fixed%use_atomic_pol
		if (verbose_mode == 1) then
			if (in_fixed%use_atomic_pol == 1) then
				write(*,FMT='(A)') ' Using atomic polarization'
			endif
			if (in_fixed%use_atomic_pol == 0) then
				write(*,FMT='(A)') ' Neglecting atomic polarization'
			endif
			if (in_fixed%use_atomic_pol == 2) then
				write(*,FMT='(A)') ' Neglecting anisotropy'
			endif
		endif		
		
! Observation angle with respect to the vertical
		call lb(12,2)
		read(12,*) in_fixed%thetad, in_fixed%chid, in_fixed%gammad
		if (verbose_mode == 1) then
			print *, 'LOS reference system'
			write(*,FMT='(4X,A,F6.2,3X,A,F6.2,3X,A,F6.2)') 'thetad : ', in_fixed%thetad, 'chid : ', in_fixed%chid, 'gammad : ', &
				in_fixed%gammad
		endif
		
! Wavelength axis
		call lb(12,2)
		read(12,*) in_fixed%omin, in_fixed%omax, in_fixed%no
		if (verbose_mode == 1) then
			print *, 'Wavelength axis'
			write(*,FMT='(4X,A,F8.3,3X,A,F8.3,3X,A,I5)') 'Min : ', in_fixed%omin, 'Max : ', in_fixed%omax, 'N. steps : ', &
				in_fixed%no
		endif
		
! Line wavelength and Doppler width
! If two components with different fields, read the two Doppler widths
		call lb(12,2)
		if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
			read(12,*) in_fixed%wl, in_params%vdopp, in_params%damping
		endif
		if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
			read(12,*) in_fixed%wl, in_params%vdopp, in_params%vdopp2, in_params%damping
		endif

! Test for the damping treatment. We can set if fixed to a number or estimated using the natural broadening
		in_fixed%damping_treatment = 0
		if (in_params%damping < 0.d0) then
			in_fixed%damping_treatment = 1
		endif

		if (verbose_mode == 1) then
			if (in_params%nslabs <= 2) then
				if (in_fixed%damping_treatment == 1) then
					write(*,FMT='(4X,A,F10.4,3X,A,F6.3,A,F6.3)') 'Wavelength [A]: ', in_fixed%wl, 'Doppler vel. [km/s] : ', &
						in_params%vdopp, 'Damping : automatic'
				else
					write(*,FMT='(4X,A,F10.4,3X,A,F6.3,A,F6.3)') 'Wavelength [A]: ', in_fixed%wl, 'Doppler vel. [km/s] : ', &
						in_params%vdopp, 'Damping : ', in_params%damping
				endif
			else
				if (in_fixed%damping_treatment == 1) then
					write(*,FMT='(4X,A,F10.4,3X,A,F6.3,A,F6.3,A,F6.3)') 'Wavelength [A]: ', in_fixed%wl, 'Doppler vel. [km/s] : ', &
						in_params%vdopp, 'Doppler vel 2. [km/s] : ', in_params%vdopp2, 'Damping : automatic'
				else
					write(*,FMT='(4X,A,F10.4,3X,A,F6.3,A,F6.3,A,F6.3)') 'Wavelength [A]: ', in_fixed%wl, 'Doppler vel. [km/s] : ', &
						in_params%vdopp, 'Doppler vel 2. [km/s] : ', in_params%vdopp2, 'Damping : ', in_params%damping
				endif
			endif
		endif
		
! Macroscopic velocity
		call lb(12,2)
		if (in_params%nslabs == 1) then
			read(12,*) in_params%vmacro
		else
			read(12,*) in_params%vmacro, in_params%vmacro2
		endif
		
! Use magneto-optical terms
		call lb(12,2)
		read(12,*) use_mag_opt_RT
		
! Use stimulated emission terms in the absorption coefficients
		call lb(12,2)
		read(12,*) use_stim_emission_RT
	
		close(12)
		
		in_fixed%total_forward_modeling = 0
		
! Take into account the observation angle (scattering angle) to calculate the height from the apparent height
		if (in_params%height < 0.d0) then
			delta = abs(90.d0 - in_fixed%thetad)
			height = (RSUN + in_params%height) / (cos(delta*PI/180.d0)) - RSUN
			write(*,FMT='(A,F7.3)') 'The height taking into account the scattering angle is : ', height
		endif
		
	
	end subroutine read_experiment
	
!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
	subroutine read_model_file(file)
	character(len=120) :: file
	integer :: i, j, nJ, n, i1, i2, ir, jmax2, jmin2, jp2, j2, kmin, kmax, k, q, loop
	real(kind=8) :: rnujjp
	character(len=128) :: fileModel(23)
				
! Generate a scratch file that will be deleted later		
		fileModel(1) = '2'
		fileModel(2) = '5'
		fileModel(3) = '1        0               '
		fileModel(4) = '                    0.00 '
		fileModel(5) = '2        2        '
      fileModel(6) = '                    0.00'
      fileModel(7) = '                   -0.987913'
      fileModel(8) = '                   -1.064340'
      fileModel(9) = '3        0       '
      fileModel(10) = '                    0.00'
      fileModel(11) = '4        2         '
      fileModel(12) = '                    0.00'
      fileModel(13) = '                   -0.270647'
      fileModel(14) = '                   -0.292616'
      fileModel(15) = ' 5        4'     
      fileModel(16) = '                    0.00'
		fileModel(17) = '                   -0.044187'
		fileModel(18) = '                   -0.046722'
		fileModel(19) = ' 4'
		fileModel(20) = '1    1    2    1.022d7    10829.0911    1.0000000    1.0000000    0.0000000'
		fileModel(21) = '2    1    4    9.478d6    3888.6046    0.2000000    1.0000000    0.0000000'
		fileModel(22) = '3    2    3    2.780d7    7065.7085    1.0000000    1.0000000    0.0000000'
		fileModel(23) = '4    2    5    7.060d7    5875.9663    1.0000000    1.0000000    0.0000000'		
		
		loop = 0

		loop = loop + 1		
		read(fileModel(loop),*) is2

		loop = loop + 1		
		read(fileModel(loop),*) n_terms		
		
		if (allocated(lsto2)) deallocate(lsto2)
		allocate(lsto2(n_terms))
		
		nrhos = 0
		
! 		open(unit=13,file='mterm.tab',action='write',status='replace')		
		
		jlimit2 = 0
		do i = 1, n_terms
			loop = loop + 1		
			read(fileModel(loop),*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			do j2 = jmin2, jmax2, 2
				loop = loop + 1		
				read(fileModel(loop),*) 
			enddo
			if (verbose_mode == 1) then
				print *, 'Level ', i
				print *, '         Jmin : ', jmin2/2.d0
				print *, '         Jmax : ', jmax2/2.d0
			endif
			jlimit2 = max(jmax2,jlimit2)
		enddo
		
		if (allocated(energy)) deallocate(energy)
		allocate(energy(n_terms,0:jlimit2))
		
		loop = 2
		
		do i = 1, n_terms
			loop = loop + 1		
			read(fileModel(loop),*) n, lsto2(i)
			
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				loop = loop + 1		
				read(fileModel(loop),*) energy(i,j2)				
				
			enddo
			
			do j2 = jmin2, jmax2, 2
				do jp2 = j2, jmax2, 2
					rnujjp = PC * (energy(i,j2)-energy(i,jp2))
					kmin = (jp2-j2) / 2
					kmax = (jp2+j2) / 2
					do k = kmin, kmax
						do q = -k, k
							if (j2 == jp2 .and. q < 0) then
							else
								do ir = 1, 2
									if (j2 == jp2 .and. q == 0 .and. ir == 2) then
									else
										nrhos = nrhos + 1
									endif
								enddo
							endif
						enddo
					enddo
				enddo
			enddo
		enddo
		
		loop = 2

		if (allocated(ntab)) deallocate(ntab)
		if (allocated(j2tab)) deallocate(j2tab)
		if (allocated(jp2tab)) deallocate(jp2tab)
		if (allocated(ktab)) deallocate(ktab)
		if (allocated(qtab)) deallocate(qtab)
		if (allocated(irtab)) deallocate(irtab)
		if (allocated(rnutab)) deallocate(rnutab)
		
		allocate(ntab(nrhos))
		allocate(j2tab(nrhos))
		allocate(jp2tab(nrhos))
		allocate(ktab(nrhos))
		allocate(qtab(nrhos))
		allocate(irtab(nrhos))
		allocate(rnutab(nrhos))
		
		nrhos = 0
		
		do i = 1, n_terms
			loop = loop + 1
			read(fileModel(loop),*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				loop = loop + 1
				read(fileModel(loop),*) energy(i,j2)
			enddo
			
			do j2 = jmin2, jmax2, 2
				do jp2 = j2, jmax2, 2
					rnujjp = PC * (energy(i,j2)-energy(i,jp2))
					kmin = (jp2-j2) / 2
					kmax = (jp2+j2) / 2
					do k = kmin, kmax
						do q = -k, k
							if (j2 == jp2 .and. q < 0) then
							else
								do ir = 1, 2
									if (j2 == jp2 .and. q == 0 .and. ir == 2) then
									else
										nrhos = nrhos + 1
										ntab(nrhos) = n
										j2tab(nrhos) = j2
										jp2tab(nrhos) = jp2
										ktab(nrhos) = k
										qtab(nrhos) = q
										irtab(nrhos) = ir
										rnutab(nrhos) = rnujjp
! 										write(13,FMT='(7(1X,I5),6X,E12.5)') nrhos, n, j2, jp2, k, q, ir, rnujjp
									endif
								enddo
							endif
						enddo
					enddo
				enddo
			enddo
		enddo

		
		! Now read transitions
		loop = loop + 1
		read(fileModel(loop),*) atom%ntran

		allocate(atom%nterml(atom%ntran))
		allocate(atom%ntermu(atom%ntran))
		allocate(atom%ae(atom%ntran))
		allocate(atom%wavelength(atom%ntran))
		allocate(atom%reduction_factor(atom%ntran))
		allocate(atom%reduction_factor_omega(atom%ntran))
		allocate(atom%j10(atom%ntran))
		
		do i = 1, atom%ntran
			loop = loop + 1
			read(fileModel(loop),*) k, atom%nterml(i), atom%ntermu(i), atom%ae(i), atom%wavelength(i), &
				atom%reduction_factor(i), atom%reduction_factor_omega(i), atom%j10(i)
		enddo
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
	end subroutine read_model_file
	

!------------------------------------------------------------
! Clean all the allocated memory
!------------------------------------------------------------	
	subroutine clean
		
		deallocate(lsto2)
		deallocate(energy)		
		deallocate(ntab)
		deallocate(ktab)
		deallocate(qtab)
		deallocate(irtab)
		deallocate(j2tab)
		deallocate(jp2tab)
		deallocate(rnutab)		
		
	end subroutine clean

!------------------------------------------------------------
! Number of columns in a line
!------------------------------------------------------------
	function ncolumns(string)
	integer :: ncolumns
	character(len=*) :: string
	character(len=1) :: ch, chprev
	integer :: i, l, blanks
		l = len(trim(string))
		blanks = 1
		chprev = string(1:1)
		do i = 2, l
			ch = string(i:i)			
			if (ch /= ' ' .and. chprev == ' ') then
				blanks = blanks + 1
			endif			
			chprev = ch
		enddo
		ncolumns = blanks - 1
	end function ncolumns
end module io
