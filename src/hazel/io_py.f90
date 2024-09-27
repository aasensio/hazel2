module io
use vars  !EDGAR: vars was commented! Check WHY. atom var structure was not being read!! 
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
		synthesis_method = 5		
		working_mode = 1
		
		
	end subroutine read_config
	


!------------------------------------------------------------
! Read any atom model file.
!------------------------------------------------------------
	subroutine read_model_file_ok(file)
	character(len=120) :: file   !EDGAR: check if it is inconsistent with the adaptative nchar length in c_init
	integer :: i, j, nJ, n, i1, i2, ir, jmax2, jmin2, jp2, j2, kmin, kmax, k, q
	real(kind=8) :: rnujjp
	!logical, save :: first_entry = .true. !EDGAR auxiliary variable

		open(unit=12,file=file,action='read',status='old')
		
		read(12,*) is2
		read(12,*) n_terms		
		
		
		if (allocated(lsto2)) deallocate(lsto2) 
		allocate(lsto2(n_terms))
		
		
		nrhos = 0
		
!		open(unit=13,file='mterm.tab',action='write',status='replace')		
		
		jlimit2 = 0
		do i = 1, n_terms
			read(12,*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			do j2 = jmin2, jmax2, 2
				read(12,*) 
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
		


		rewind(12)
		call lb(12,2)
		file_pointer = 2
		
		do i = 1, n_terms
			read(12,*) n, lsto2(i)
			file_pointer = file_pointer + 1
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				read(12,*) energy(i,j2)
				file_pointer = file_pointer + 1
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
		
		rewind(12)
		call lb(12,2)

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
			read(12,*) n, lsto2(i)
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
			do j2 = jmin2, jmax2, 2
				read(12,*) energy(i,j2)
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
!										write(13,FMT='(7(1X,I5),6X,E12.5)') nrhos, n, j2, jp2, k, q, ir, rnujjp
									endif
								enddo
							endif
						enddo
					enddo
				enddo
			enddo
		enddo

! Now read transitions

!EDGAR: the next block cannot be done because pointers inside atom are no allocatable, are just free
!pointers and what you have to allocate are the variables where they will point
! 		if (allocated(atom%nterml)) deallocate(atom%nterml)
! 		if (allocated(atom%ntermu)) deallocate(atom%ntermu)
! 		if (allocated(atom%ae)) deallocate(atom%ae)
! 		if (associated(atom%wavelength)) deallocate(atom%wavelength)
! 		if (allocated(atom%reduction_factor)) deallocate(atom%reduction_factor)
! 		if (allocated(atom%reduction_factor_omega)) deallocate(atom%reduction_factor_omega)
! 		if (allocated(atom%j10)) deallocate(atom%j10)


		read(12,*) atom%ntran

		allocate(atom%nterml(atom%ntran))
		allocate(atom%ntermu(atom%ntran))
		allocate(atom%ae(atom%ntran))
		!EDGAR . nullify the pointer is not really needed because we inmediately allocate 
        nullify(atom%wavelength)    
       	allocate(atom%wavelength(atom%ntran))
		allocate(atom%reduction_factor(atom%ntran))
		allocate(atom%reduction_factor_omega(atom%ntran))

		nullify(atom%j10)       !EDGAR ...same as above 
        allocate(atom%j10(atom%ntran))   
		!if (.not. associated(atom%j10)) allocate(atom%j10(atom%ntran))
		
		!here j10 is readed from file, but it will be overwritten later by the values entered from 
		!python routine
		do i = 1, atom%ntran
			read(12,*) k, atom%nterml(i), atom%ntermu(i), atom%ae(i), atom%wavelength(i), &
				atom%reduction_factor(i), atom%reduction_factor_omega(i), atom%j10(i)
			!print*,'atom%reduction_factor',atom%reduction_factor(i) 
			!print*,'j10 from io_py.f90:',atom%j10(i)
		enddo

		close(12)
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		return  
	end subroutine read_model_file_ok
	
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
		is2 = 2

		loop = loop + 1		
		n_terms = 5
		
		if (allocated(lsto2)) deallocate(lsto2)
		allocate(lsto2(n_terms))

		lsto2(1) = 0
		lsto2(2) = 2
		lsto2(3) = 0
		lsto2(4) = 2
		lsto2(5) = 4
		
		nrhos = 0
		
! 		open(unit=13,file='mterm.tab',action='write',status='replace')		
		
		jlimit2 = 0
		do i = 1, n_terms
			loop = loop + 1		
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			! do j2 = jmin2, jmax2, 2
			! 	loop = loop + 1		
			! 	print *, i, j2, fileModel(loop)
				
			! enddo
			if (verbose_mode == 1) then
				print *, 'Level ', i
				print *, '         Jmin : ', jmin2/2.d0
				print *, '         Jmax : ', jmax2/2.d0
			endif
			jlimit2 = max(jmax2,jlimit2)
		enddo
		
		if (allocated(energy)) deallocate(energy)
		allocate(energy(n_terms,0:jlimit2))

		energy(1,2) = 0.0
        energy(2,0) = 0.0
        energy(2,2) = -0.987913
        energy(2,4) = -1.064340
        energy(3,2) = 0.0
        energy(4,0) = 0.0
        energy(4,2) = -0.270647
        energy(4,4) = -0.292616
        energy(5,2) = 0.0
        energy(5,4) = -0.044187
        energy(5,6) = -0.046722
		
		loop = 2
		
		do i = 1, n_terms
			
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
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
			n = i
			
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
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
		! loop = loop + 1
		! read(fileModel(loop),*) atom%ntran

		atom%ntran = 4

		allocate(atom%nterml(atom%ntran))
		allocate(atom%ntermu(atom%ntran))
		allocate(atom%ae(atom%ntran))
		allocate(atom%wavelength(atom%ntran))
		allocate(atom%reduction_factor(atom%ntran))
		allocate(atom%reduction_factor_omega(atom%ntran))
		allocate(atom%j10(atom%ntran))

		atom%nterml(1) = 1
		atom%ntermu(1) = 2
		atom%ae(1) = 1.022d7
		atom%wavelength(1) = 10829.0911d0
		atom%reduction_factor(1) = 1.d0
		atom%reduction_factor_omega(1) = 1.d0
		atom%j10(1) = 0.d0

		atom%nterml(2) = 1
		atom%ntermu(2) = 4
		atom%ae(2) = 9.478d6
		atom%wavelength(2) = 3888.6046d0
		atom%reduction_factor(2) = 0.2d0
		atom%reduction_factor_omega(2) = 1.0d0
		atom%j10(2) = 0.0d0

		atom%nterml(3) = 2
		atom%ntermu(3) = 3
		atom%ae(3) = 2.780d7
		atom%wavelength(3) = 7065.7085d0
		atom%reduction_factor(3) = 1.0d0
		atom%reduction_factor_omega(3) = 1.0d0
		atom%j10(3) = 0.0d0

		atom%nterml(4) = 2
		atom%ntermu(4) = 5
		atom%ae(4) = 7.060d7
		atom%wavelength(4) = 5875.9663d0
		atom%reduction_factor(4) = 1.0d0
		atom%reduction_factor_omega(4) = 1.0d0
		atom%j10(4) = 0.0d0


		! do i = 1, atom%ntran
		! 	loop = loop + 1
		! 	read(fileModel(loop),*) k, atom%nterml(i), atom%ntermu(i), atom%ae(i), atom%wavelength(i), &
		! 		atom%reduction_factor(i), atom%reduction_factor_omega(i), atom%j10(i)
		! enddo
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
	end subroutine read_model_file
	

!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
	subroutine read_model_file_sodium(file)
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
		is2 = 3

		loop = loop + 1		
		n_terms = 3
		
		if (allocated(lsto2)) deallocate(lsto2)
		allocate(lsto2(n_terms))

		lsto2(1) = 1
		lsto2(2) = 1
		lsto2(3) = 3
		
		nrhos = 0
		
! 		open(unit=13,file='mterm.tab',action='write',status='replace')		
		
		jlimit2 = 0
		do i = 1, n_terms
			loop = loop + 1		
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			! do j2 = jmin2, jmax2, 2
			! 	loop = loop + 1		
			! 	print *, i, j2, fileModel(loop)
				
			! enddo
			if (verbose_mode == 1) then
				print *, 'Level ', i
				print *, '         Jmin : ', jmin2/2.d0
				print *, '         Jmax : ', jmax2/2.d0
			endif
			jlimit2 = max(jmax2,jlimit2)
		enddo
		
		if (allocated(energy)) deallocate(energy)
		allocate(energy(n_terms,0:jlimit2))

		energy(1,2) = -0.036906
		energy(1,4) = 0.0221437
        energy(2,2) = -0.003929
        energy(2,4) = 0.002357		
        energy(3,0) = -0.002325
        energy(3,2) = -0.001705
        energy(3,4) = -0.000465
        energy(3,6) = 0.001395
		
		loop = 2
		
		do i = 1, n_terms
			
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
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
			n = i
			
			jmin2 = abs(is2-lsto2(i))
			jmax2 = is2 + lsto2(i)
			
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
		! loop = loop + 1
		! read(fileModel(loop),*) atom%ntran

		atom%ntran = 2

		allocate(atom%nterml(atom%ntran))
		allocate(atom%ntermu(atom%ntran))
		allocate(atom%ae(atom%ntran))
		allocate(atom%wavelength(atom%ntran))
		allocate(atom%reduction_factor(atom%ntran))
		allocate(atom%reduction_factor_omega(atom%ntran))
		allocate(atom%j10(atom%ntran))

		atom%nterml(1) = 1
		atom%ntermu(1) = 2
		atom%ae(1) = 6.14d7
		atom%wavelength(1) = 5895.924d0
		atom%reduction_factor(1) = 1.d0
		atom%reduction_factor_omega(1) = 1.d0
		atom%j10(1) = 0.d0

		atom%nterml(2) = 1
		atom%ntermu(2) = 3
		atom%ae(2) = 6.16d7
		atom%wavelength(2) = 5889.950d0
		atom%reduction_factor(2) = 1.d0
		atom%reduction_factor_omega(2) = 1.0d0
		atom%j10(2) = 0.0d0
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
		if (verbose_mode == 1) then
			print *, 'Number of unknowns : ', nrhos
		endif
		
	end subroutine read_model_file_sodium


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
