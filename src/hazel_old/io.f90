module io
! use vars
use maths
#if defined(mpi)
use netcdf
#endif
implicit none
contains

!------------------------------------------------------------
! Read the main configuration file
!------------------------------------------------------------
    subroutine read_config
    character(len=100) :: config_file, buffer
    integer :: nargs, ios, iargc

! Verify if a configuration file is passed as an argument
! If not, use the default 'config' file

! Call the code with
! phazel config_inversion.dat 1 5000

        nargs = iargc()

        if (nargs == 0) then
            config_file = 'config_inversion.dat'
            starting_pixel = 1
            final_pixel = -99
        endif
        
        if (nargs == 1) then
            call getarg(1,config_file)
            starting_pixel = 1
            final_pixel = -99
        endif

        if (nargs == 2) then
            call getarg(1,config_file)
            call getarg(2,buffer)
            read(buffer,*) starting_pixel
            print *, 'Starting from pixel ', starting_pixel
            final_pixel = -99
        endif

        if (nargs == 3) then
            call getarg(1,config_file)
            call getarg(2,buffer)
            read(buffer,*) starting_pixel
            call getarg(3,buffer)
            read(buffer,*) final_pixel
            print *, 'Starting from pixel ', starting_pixel
            print *, 'Ending at pixel ', final_pixel
        endif
    
        open(unit=12,file=config_file,action='read',status='old',IOSTAT=ios)
        if (ios /= 0) then
            write(*,*) 'Configuration file not present. Verify that a configuration file is present'
            stop
        endif
        call lb(12,5)
        read(12,*) input_model_file
        call lb(12,2)
        read(12,*) input_experiment
        call lb(12,2)
        read(12,*) direct_ranges
        call lb(12,2)
        read(12,*) output_rho_vertical_upper
        call lb(12,2)
        read(12,*) output_rho_vertical_lower
        call lb(12,2)
        read(12,*) output_rho_magnetic_upper
        call lb(12,2)
        read(12,*) output_rho_magnetic_lower
        call lb(12,2)
        read(12,*) output_rtcoef
        call lb(12,2)
        read(12,*) output_rtcoef_zeeman
        call lb(12,2)
        read(12,*) input_observed_profiles
        call lb(12,2)
        read(12,*) output_inverted_profiles
        call lb(12,2)
        read(12,*) output_final_parameters
        call lb(12,2)
        read(12,*) output_error_parameters
        call lb(12,2)
        read(12,*) input_inverted_parameters
        call lb(12,2)
        read(12,*) verbose_mode
        call lb(12,2)
        read(12,*) linear_solver
        call lb(12,2)
        read(12,*) synthesis_mode
        call lb(12,2)
        read(12,*) working_mode     
        close(12)
        
    end subroutine read_config
    

!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
    subroutine read_experiment(file,in_params,in_fixed)
    character(len=120) :: file
    type(variable_parameters) :: in_params
    type(fixed_parameters) :: in_fixed
    real(kind=8) :: delta, height
    integer :: j, i
    
        nparam = 0
        
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
            nparam = nparam + 6
        else
            read(12,*) in_params%bgauss, in_params%thetabd, in_params%chibd
            nparam = nparam + 3
        endif
        if (imag == 1 .and. verbose_mode == 1) then
            print *, 'Field properties'
            write(*,FMT='(4X,A,F8.3,3X,A,F6.2,3X,A,F6.2)') 'B : ', in_params%bgauss, 'thetab : ', in_params%thetabd, 'chib : ',&
                in_params%chibd
        endif       
        
! Values of the height for the calculation of the radiation anisotropy
        call lb(12,2)
        read(12,*) in_params%height

        nparam = nparam + 1
                
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
            nparam = nparam + 1
        endif

        if (in_params%nslabs == 2 .or. in_params%nslabs == 3) then
            read(12,*) in_params%dtau, in_params%dtau2
            nparam = nparam + 2
        endif

        if (in_params%nslabs == -2) then
            read(12,*) in_params%dtau, in_params%dtau2, in_params%ff
            nparam = nparam + 3
        endif
        
! Value of the gradient of the source function (ME)
        call lb(12,2)
        if (in_params%nslabs == 1) then
            read(12,*) in_params%beta
            nparam = nparam + 1
        else
            read(12,*) in_params%beta, in_params%beta2
            nparam = nparam + 2
        endif

        
! Incident Stokes parameters
        call lb(12,2)
        read(12,*) in_fixed%Stokes_incident_mode        
        if (index(in_fixed%Stokes_incident_mode, 'single') /= 0) then
            read(12,*) (in_fixed%Stokes_incident(j),j=0,3)
        endif
        if (index(in_fixed%Stokes_incident_mode, 'file') /= 0) then
            read(12,*) in_fixed%Stokes_incident_file
            open(unit=18, file=in_fixed%Stokes_incident_file, action='read', status='old')
            read(18,*) in_fixed%stokes_boundary_len
            allocate(in_fixed%stokes_boundary(0:3,in_fixed%stokes_boundary_len))
            do i = 1, in_fixed%stokes_boundary_len
                read(18,*) (in_fixed%stokes_boundary(j,i),j=0,3)
            enddo
            close(18)
        endif
                
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
        call lb(12,2)
        if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
            read(12,*) in_fixed%wl, in_params%vdopp, in_params%damping
            nparam = nparam + 2
        endif
        if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            read(12,*) in_fixed%wl, in_params%vdopp, in_params%vdopp2, in_params%damping
            nparam = nparam + 3
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
            nparam = nparam + 1
        else
            read(12,*) in_params%vmacro, in_params%vmacro2
            nparam = nparam + 2
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
    integer :: i, j, nJ, n, i1, i2, ir, jmax2, jmin2, jp2, j2, kmin, kmax, k, q
    real(kind=8) :: rnujjp

        open(unit=12,file=file,action='read',status='old')
        
        read(12,*) is2
        read(12,*) n_terms      
        
        allocate(lsto2(n_terms))
        
        nrhos = 0
        
!       open(unit=13,file='mterm.tab',action='write',status='replace')      
        
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
!                                       write(13,FMT='(7(1X,I5),6X,E12.5)') nrhos, n, j2, jp2, k, q, ir, rnujjp
                                    endif
                                enddo
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

! Now read transitions
        read(12,*) atom%ntran

        allocate(atom%nterml(atom%ntran))
        allocate(atom%ntermu(atom%ntran))
        allocate(atom%ae(atom%ntran))
        allocate(atom%wavelength(atom%ntran))
        allocate(atom%reduction_factor(atom%ntran))
        allocate(atom%reduction_factor_omega(atom%ntran))
        allocate(atom%j10(atom%ntran))
        
        do i = 1, atom%ntran
            read(12,*) k, atom%nterml(i), atom%ntermu(i), atom%ae(i), atom%wavelength(i), &
                atom%reduction_factor(i), atom%reduction_factor_omega(i), atom%j10(i)
        enddo

        close(12)
        
        if (verbose_mode == 1) then
            print *, 'Number of unknowns : ', nrhos
        endif
        
    end subroutine read_model_file
    

!------------------------------------------------------------
! Read the DIRECT configuration file
!------------------------------------------------------------
    subroutine read_direct_conf(direct_ranges, in_fixed)
    character(len=120) :: direct_ranges
    type(fixed_parameters) :: in_fixed
    integer :: i
        
        allocate(in_fixed%upper_direct(18))
        allocate(in_fixed%lower_direct(18))

        open(unit=24,file=direct_ranges,action='read',status='old')
        call lb(24,5)

        read(24,*) in_fixed%DIRmaxf
        if (in_fixed%DIRmaxf < 0) in_fixed%DIRmaxf = 1000
!       print *, 'Maximum number of function evaluations : ', in_fixed%DIRmaxf
        
        call lb(24,2)
        read(24,*) in_fixed%volper
        if (in_fixed%volper < 0) in_fixed%volper = 1.d-20
!       print *, 'Stopping when the hypervolume with respect to the original is less than ', in_fixed%volper
        
        do i = 1, 18
            call lb(24,2)
            read(24,*) in_fixed%lower_direct(i), in_fixed%upper_direct(i)
        enddo
        close(24)
        
    end subroutine read_direct_conf

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
! Set boundary condition
!------------------------------------------------------------
    subroutine set_boundary_condition(in_fixed, in_inversion)
    type(fixed_parameters) :: in_fixed
    type(type_inversion) :: in_inversion
    integer :: j

        if (index(in_fixed%Stokes_incident_mode, 'single') /= 0) then
            if (.not. associated(in_fixed%stokes_boundary)) allocate(in_fixed%stokes_boundary(0:3,in_fixed%no))
            do j = 0, 3
                in_fixed%stokes_boundary(j,:) = in_fixed%Stokes_incident(j)
            enddo
        endif
        if (index(in_fixed%Stokes_incident_mode, 'file') /= 0) then
            if (in_fixed%stokes_boundary_len /= in_fixed%no) then
                write(*,*) 'Wavelength axis in boundary condition file is different from synthesis wavelength axis'
                stop
            endif
        endif

    end subroutine set_boundary_condition

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

!------------------------------------------------------------
! Check if the NetCDF action ended correctly
!------------------------------------------------------------
    subroutine check(status)
   integer, intent ( in) :: status

#if defined(mpi)
        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
        stop
        end if
#endif
    endsubroutine check
  
!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
    subroutine preread_observations(file_obs,in_fixed,in_observation)
    character(len=120) :: file_obs
    character(len=320) :: column, name_var
    type(type_observation) :: in_observation
    type(fixed_parameters) :: in_fixed
    integer :: i, j, format_file, res, number_columns
    real(kind=8) :: temp
    integer :: ncid, npixel, ncol, nlambda, test(2)

! Test if this is a single file (ending in .prof) or a NetCDF file with a map
        if (index(file_obs, '.prof') /= 0) then
            
            open(unit=12,file=file_obs,action='read',status='old')
            read(12,*) in_observation%n, in_observation%normalization

! Now detect the file type (4 columns for a simple estimation of sigma from the amplitude)
! and 8 columns for explicitly giving the value of sigma for each wavelength)
            read(12,FMT='(A)') column
            number_columns = ncolumns(column)
            close(12)
            if (number_columns == 5) then
                in_observation%observation_format = 0
            endif
            if (number_columns == 9) then
                in_observation%observation_format = 1
            endif
            close(12)

! Only 1 profile is considered
            in_observation%npixel = 1
        endif
        
        if (index(file_obs, '.nc') /= 0) then
#if defined(mpi)
            in_observation%observation_format = 2
            
! Open the file in read-only access
            call check( nf90_open(file_obs, NF90_NOWRITE, in_observation%obs_id) )

! Get the varid of the data variable, based on its name.
            call check( nf90_inq_dimid(in_observation%obs_id, "npixel", in_observation%pix_id) )
            call check( nf90_inq_dimid(in_observation%obs_id, "ncolumns", in_observation%col_id) )
            call check( nf90_inq_dimid(in_observation%obs_id, "nstokes_par", in_observation%nstokespar_id) )
            call check( nf90_inq_dimid(in_observation%obs_id, "nlambda", in_observation%nlambda_id) )
            
            call check( nf90_inq_varid(in_observation%obs_id, "normalization", in_observation%normalization_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "map", in_observation%map_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "lambda", in_observation%lambda_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "boundary", in_observation%boundary_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "height", in_observation%height_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "obs_theta", in_observation%obstheta_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "obs_gamma", in_observation%obsgamma_id) )
            call check( nf90_inq_varid(in_observation%obs_id, "pars_initial", in_observation%parsInit_id) )
            
            
            call check( nf90_inquire_dimension(in_observation%obs_id, in_observation%pix_id, name_var, npixel) )
            call check( nf90_inquire_dimension(in_observation%obs_id, in_observation%col_id, name_var, ncol) )
            call check( nf90_inquire_dimension(in_observation%obs_id, in_observation%nlambda_id, name_var, nlambda) )

            call check( nf90_get_var(in_observation%obs_id, in_observation%normalization_id, temp) )

            if (temp == 0) then
                in_observation%normalization = 'cont'
            endif
            if (temp == 1) then
                in_observation%normalization = 'peak'
            endif
                        
            if (verbose_mode == 1) then
                write(*,FMT='(A,I6)') 'Number of pixels : ', npixel
                write(*,FMT='(A,I4)') 'Number of wavelengths : ', nlambda
                write(*,FMT='(A,A4)') 'Normalization : ', in_observation%normalization
            endif

! Now that we have opened the file with the observations, open two files with the
! results and the inverted profiles

! If start from the beginning
            if (starting_pixel == 1) then
                call check( nf90_create(output_inverted_profiles, NF90_CLOBBER, in_fixed%syn_id) )
                call check( nf90_create(output_final_parameters, NF90_CLOBBER, in_fixed%par_id) )
                call check( nf90_create(output_error_parameters, NF90_CLOBBER, in_fixed%error_id) )
                
! Set dimensions and set variables for the synthetic profiles
                call check( nf90_def_dim(in_fixed%syn_id, "npixel", npixel, in_fixed%pix_syn_id) )
                call check( nf90_def_dim(in_fixed%syn_id, "ncolumns", 4, in_fixed%col_syn_id) )
                call check( nf90_def_dim(in_fixed%syn_id, "nlambda", nlambda, in_fixed%nlambda_syn_id) )
                
                call check( nf90_def_var(in_fixed%syn_id, "lambda", NF90_DOUBLE, &
                    (/ in_fixed%nlambda_syn_id /), in_fixed%lambda_syn_id) )
                call check( nf90_def_var(in_fixed%syn_id, "map", NF90_DOUBLE, &
                    (/ in_fixed%col_syn_id, in_fixed%nlambda_syn_id, in_fixed%pix_syn_id /), in_fixed%map_syn_id) )
                call check( nf90_enddef(in_fixed%syn_id) )
                
! Set dimensions and set variables for the inverted values
                call check( nf90_def_dim(in_fixed%par_id, "npixel", npixel, in_fixed%pix_par_id) )
                call check( nf90_def_dim(in_fixed%par_id, "ncolumns", nparam, in_fixed%col_par_id) )
                call check( nf90_def_var(in_fixed%par_id, "map", NF90_DOUBLE, (/ in_fixed%col_par_id, in_fixed%pix_par_id/), in_fixed%map_par_id) )
                call check( nf90_enddef(in_fixed%par_id) )
                
! Set dimensions and set variables for the inverted values
                call check( nf90_def_dim(in_fixed%error_id, "npixel", npixel, in_fixed%pix_error_id) )
                call check( nf90_def_dim(in_fixed%error_id, "ncolumns", nparam, in_fixed%col_error_id) )
                call check( nf90_def_var(in_fixed%error_id, "map", NF90_DOUBLE, (/ in_fixed%col_error_id, in_fixed%pix_error_id/), in_fixed%map_error_id) )
                call check( nf90_enddef(in_fixed%error_id) )

            endif

            if (starting_pixel /= 1) then
! Open file with results
                call check( nf90_open(output_inverted_profiles, NF90_WRITE, in_fixed%syn_id) )
                call check( nf90_open(output_final_parameters, NF90_WRITE, in_fixed%par_id) )
                call check( nf90_open(output_error_parameters, NF90_WRITE, in_fixed%error_id) )

! Inquire dimensions
                call check( nf90_inq_dimid(in_fixed%syn_id, "npixel", in_fixed%pix_syn_id) )
                call check( nf90_inq_dimid(in_fixed%syn_id, "ncolumns", in_fixed%col_syn_id) )
                call check( nf90_inq_dimid(in_fixed%syn_id, "nlambda", in_fixed%nlambda_syn_id) )
                call check( nf90_inq_dimid(in_fixed%par_id, "npixel", in_fixed%pix_par_id) )
                call check( nf90_inq_dimid(in_fixed%par_id, "ncolumns", in_fixed%col_par_id) )
                call check( nf90_inq_dimid(in_fixed%error_id, "npixel", in_fixed%pix_error_id) )
                call check( nf90_inq_dimid(in_fixed%error_id, "ncolumns", in_fixed%col_error_id) )

! Define variables
                call check( nf90_inq_varid(in_fixed%syn_id, "lambda", in_fixed%lambda_syn_id) )
                call check( nf90_inq_varid(in_fixed%syn_id, "map", in_fixed%map_syn_id) )
                call check( nf90_inq_varid(in_fixed%par_id, "map", in_fixed%map_par_id) )
                call check( nf90_inq_varid(in_fixed%error_id, "map", in_fixed%map_error_id) )

! Get some dimensions
                call check( nf90_inquire_dimension(in_fixed%syn_id, in_fixed%nlambda_syn_id, name_var, nlambda) )
                call check( nf90_inquire_dimension(in_fixed%syn_id, in_fixed%pix_syn_id, name_var, npixel) )
                
            endif

! Number of wavelengths
            in_observation%n = nlambda
            in_observation%npixel = npixel
#else
    print *, 'NetCDF not supported in serial version'   
#endif
        endif

! If no final pixel selected, run until the end
        if (final_pixel == -99) then
            final_pixel = in_observation%npixel
        endif

! In the inversion mode, the wavelength axis has to be the same as that introduced in the observation
        if (working_mode == 1) then
            in_fixed%no = in_observation%n
        endif

    end subroutine preread_observations

!------------------------------------------------------------
! Read the file with the data for the experiment
!------------------------------------------------------------
    subroutine read_observation(file_obs,in_fixed,in_observation,in_params,pixel,status)
    character(len=120) :: file_obs
    character(len=320) :: column
    type(fixed_parameters) :: in_fixed
    type(variable_parameters) :: in_params
    type(type_observation) :: in_observation
    integer :: i, j, format_file, pixel, status
    real(kind=8) :: SI_max, SQ_max, SU_max, SV_max, Noise, MinAmpli
    real(kind=8), allocatable :: values(:,:), values_vec(:)
    integer, allocatable :: start(:), count(:)
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    integer :: valuesTime(8), nParamRead

! If we are doing one-profile inversion, just read data
        if (in_observation%observation_format == 0 .or. in_observation%observation_format == 1) then

! Re-read the first line of the file
            open(unit=12,file=file_obs,action='read',status='old')
            read(12,*) in_observation%n

            if (in_observation%observation_format == 0) then
                write(*,*) 'Using noise levels to equalize the amplitudes of all Stokes parameters'
                do i = 1, in_observation%n
                    read(12,*) in_observation%wl(i), (in_observation%stokes(j,i),j=0,3)
                enddo
                SI_max = maxval(abs(in_observation%stokes(0,:)))
                SQ_max = maxval(abs(in_observation%stokes(1,:)))
                SU_max = maxval(abs(in_observation%stokes(2,:)))
                SV_max = maxval(abs(in_observation%stokes(3,:)))
            
                in_observation%sigma(0,:) = SI_max
                in_observation%sigma(1,:) = SQ_max
                in_observation%sigma(2,:) = SU_max
                in_observation%sigma(3,:) = SV_max
            
                Noise = 1.d-4
                if (SI_max < 1.d-20) SI_max = 1.d0
                if (SQ_max < Noise*20) SQ_max = SI_max / 50.d0
                if (SU_max < Noise*20) SU_max = SI_max / 50.d0
                if (SV_max < Noise*20) SV_max = SI_max / 20.d0
            
                MinAmpli = min(SI_max,SQ_max,SU_max,SV_max)

                in_observation%sigma = in_observation%sigma / sqrt(1.d3)
            endif

            if (in_observation%observation_format == 1) then
                    write(*,*) 'Reading noise levels from file'
                do i = 1, observation%n
                    read(12,*) in_observation%wl(i), (in_observation%stokes(j,i),j=0,3), (in_observation%sigma(j,i),j=0,3)
                enddo
            endif
        
            close(12)

            status = 0
            
        endif

! If we are inverting a map, read the appropriate point and put it onto the arrays
        if (in_observation%observation_format == 2) then            
#if defined(mpi)
            if (pixel <= final_pixel) then
            
                call date_and_time(date, time, zone, valuesTime)
                write(*,FMT='(A10,A,I6,A,I6)') time, ' * Master - Reading pixel ', pixel, ' from ', final_pixel
                
                allocate(values(8,in_observation%n))
                
                allocate(start(3))
                allocate(count(3))
                start = (/ 1, 1, pixel /)
                count = (/ 8, in_observation%n, 1 /)
                
                call date_and_time(date, time, zone, valuesTime)

                
! Read data for the appropriate pixel               
                call check( nf90_get_var(in_observation%obs_id, in_observation%lambda_id, in_observation%wl) )
                call check( nf90_get_var(in_observation%obs_id, in_observation%map_id, values, start=start, count=count) )
                                                
                call date_and_time(date, time, zone, valuesTime)

! Put the profile on the arrays             
                in_observation%stokes(0:3,:) = values(1:4,:)
                in_observation%sigma(0:3,:) = values(5:8,:)

                deallocate(values, start, count)
                
                call date_and_time(date, time, zone, valuesTime)
                                
! Read incident Stokes parameter            
                call check( nf90_get_var(in_observation%obs_id, in_observation%boundary_id, in_fixed%Stokes_incident,&
                    start=(/ 1, pixel /), count=(/ 4, 1 /)) )

                call date_and_time(date, time, zone, valuesTime)                
                                    
! Read height
                allocate(values_vec(1))             
                call check( nf90_get_var(in_observation%obs_id, in_observation%height_id, values_vec,&
                    start=(/ pixel /), count=(/ 1 /) ) )
                    
                call date_and_time(date, time, zone, valuesTime)
                
                in_params%height = values_vec(1)
                                
! Read observation theta angle              
                call check( nf90_get_var(in_observation%obs_id, in_observation%obstheta_id, values_vec,&
                    start=(/ pixel /), count=(/ 1 /) ) )
                    
                call date_and_time(date, time, zone, valuesTime)
                
                in_fixed%thetad = values_vec(1)

! Read observation reference gamma angle                
                call check( nf90_get_var(in_observation%obs_id, in_observation%obsgamma_id, values_vec,&
                    start=(/ pixel /), count=(/ 1 /) ) )
                    
                call date_and_time(date, time, zone, valuesTime)

                in_fixed%gammad = values_vec(1)

! Read the initial values for the parameters. Those that are inverted will be modified (starting from these values if using LM)
! but those that are kept fixed
                deallocate(values_vec)

! 1-component (vector of size 8): B, thetaB, chiB, tau, vdop, a, vmac, beta
! 2-component 1+1 with same field (vector of size 11): B, thetaB, chiB, tau1, tau2, vdop, a, vmac1, vmac2, beta, beta2
! 2-component 1+1 with different field (vector of size 15): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, beta, beta2
! 2-component 2 with different field with ff (vector of size 16): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, ff, beta1, beta2
                select case(in_params%nslabs)
                    case(1)
                        allocate(values_vec(8))
                        nParamRead = 8
                    case(2)
                        allocate(values_vec(11))
                        nParamRead = 11
                    case(3)
                        allocate(values_vec(15))
                        nParamRead = 15
                    case(-2)
                        allocate(values_vec(16))
                        nParamRead = 16
                end select
                
                call check( nf90_get_var(in_observation%obs_id, in_observation%parsInit_id, values_vec,&
                    start=(/ 1, pixel /), count=(/ nParamRead, 1 /)) )

                call date_and_time(date, time, zone, valuesTime)
                
! Set all initial values, except for the height, which is already set before
! Single component case
                if (in_params%nslabs == 1) then
                    in_params%bgauss = values_vec(1)
                    in_params%thetabd = values_vec(2)
                    in_params%chibd = values_vec(3)
                    in_params%dtau = values_vec(4)
                    in_params%vdopp = values_vec(5)
                    in_params%damping = values_vec(6)
                    in_params%vmacro = values_vec(7)
                    in_params%beta = values_vec(8)
                endif

! Two components one after the other with the same field
                if (in_params%nslabs == 2) then
                    in_params%bgauss = values_vec(1)
                    in_params%thetabd = values_vec(2)
                    in_params%chibd = values_vec(3)                 
                    in_params%dtau = values_vec(4)
                    in_params%dtau2 = values_vec(5)
                    in_params%vdopp = values_vec(6)
                    in_params%damping = values_vec(7)
                    in_params%vmacro = values_vec(8)
                    in_params%vmacro2 = values_vec(9)
                    in_params%beta = values_vec(10)
                    in_params%beta2 = values_vec(11)
                endif

! Two components one after the other with different fields
                if (in_params%nslabs == 3) then
                    in_params%bgauss = values_vec(1)
                    in_params%thetabd = values_vec(2)
                    in_params%chibd = values_vec(3)
                    in_params%bgauss2 = values_vec(4)
                    in_params%thetabd2 = values_vec(5)
                    in_params%chibd2 = values_vec(6)
                    in_params%dtau = values_vec(7)
                    in_params%dtau2 = values_vec(8)
                    in_params%vdopp = values_vec(9)
                    in_params%vdopp2 = values_vec(10)
                    in_params%damping = values_vec(11)
                    in_params%vmacro = values_vec(12)
                    in_params%vmacro2 = values_vec(13)
                    in_params%beta = values_vec(14)
                    in_params%beta2 = values_vec(15)
                endif

! Two components with filling factor
                if (in_params%nslabs == -2) then
                    in_params%bgauss = values_vec(1)
                    in_params%thetabd = values_vec(2)
                    in_params%chibd = values_vec(3)
                    in_params%bgauss2 = values_vec(4)
                    in_params%thetabd2 = values_vec(5)
                    in_params%chibd2 = values_vec(6)                    
                    in_params%dtau = values_vec(7)
                    in_params%dtau2 = values_vec(8)
                    in_params%vdopp = values_vec(9)
                    in_params%vdopp2 = values_vec(10)
                    in_params%damping = values_vec(11)
                    in_params%vmacro = values_vec(12)
                    in_params%vmacro2 = values_vec(13)
                    in_params%ff = values_vec(14)
                    in_params%beta = values_vec(15)
                    in_params%beta2 = values_vec(16)
                endif
                                
                deallocate(values_vec)
                
!               call date_and_time(date, time, zone, valuesTime)
!               write(*,FMT='(A10,A,I6,A,I6)') time, ' - Read pixel ', pixel, ' from ', final_pixel
                
                status = 1
            else
                write(*,*) 'I have nothing to do. Killing slave...'
                status = 0
            endif
#else
    print *, 'NetCDF input not supported in serial version'
#endif
        endif
                        
    end subroutine read_observation
    
!------------------------------------------------------------
! Write results in multipixel inversion
!------------------------------------------------------------
    subroutine write_results(in_fixed,in_observation,in_inversion,in_params,error,pixel)
    type(fixed_parameters) :: in_fixed
    type(type_observation) :: in_observation
    type(type_inversion) :: in_inversion
    type(variable_parameters) :: in_params, error
    integer :: pixel, status
    real(kind=8), allocatable :: values(:,:), values_vec(:)
    integer, allocatable :: start(:), count(:)

#if defined(mpi)
        if (in_observation%observation_format == 2) then
            allocate(values(4,in_observation%n))

! Save wavelength axis
            call check( nf90_put_var(in_fixed%syn_id, in_fixed%lambda_syn_id, in_observation%wl) )

! Get profiles
            values = in_inversion%stokes_unperturbed(0:3,:)                         
            
            allocate(start(3))
            allocate(count(3))
            start = (/ 1, 1, pixel /)
            count = (/ 4, in_observation%n, 1 /)
            call check( nf90_put_var(in_fixed%syn_id, in_fixed%map_syn_id, values, start=start, count=count) )

            deallocate(values, start, count)
        endif

        allocate(values_vec(nparam))

! Set the array according to the number of components
        if (nparam == 9) then
            values_vec = (/ in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%height, in_params%dtau, in_params%vdopp, &
                in_params%damping, in_params%vmacro, in_params%beta /)
        endif

        if (nparam == 12) then
            values_vec = (/ in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%height, in_params%dtau, in_params%dtau2, &
                in_params%vdopp, in_params%damping, in_params%vmacro, in_params%vmacro2, in_params%beta, in_params%beta2 /)
        endif

        if (nparam == 16) then
            values_vec = (/ in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%bgauss2, in_params%thetabd2, in_params%chibd2,&
                in_params%height, in_params%dtau, in_params%dtau2, in_params%vdopp, in_params%vdopp2, in_params%damping, &
                in_params%vmacro, in_params%vmacro2, in_params%beta, in_params%beta2 /)
        endif

        if (nparam == 17) then
            values_vec = (/ in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%bgauss2, in_params%thetabd2, in_params%chibd2,&
                in_params%height, in_params%dtau, in_params%dtau2, in_params%vdopp, in_params%vdopp2, in_params%damping, &
                in_params%vmacro, in_params%vmacro2, in_params%beta, in_params%ff, in_params%beta2 /)
        endif
            
        allocate(start(2))
        allocate(count(2))
        
        start = (/ 1, pixel /)
        count = (/ nparam, 1 /)
        call check( nf90_put_var(in_fixed%par_id, in_fixed%map_par_id, values_vec, start=start, count=count) )
        
! Set the error array according to the number of components
        if (nparam == 9) then
            values_vec = (/ error%bgauss, error%thetabd, error%chibd, error%height, error%dtau, error%vdopp, &
                error%damping, error%vmacro, error%beta /)
        endif

        if (nparam == 12) then
            values_vec = (/ error%bgauss, error%thetabd, error%chibd, error%height, error%dtau, error%dtau2, &
                error%vdopp, error%damping, error%vmacro, error%vmacro2, error%beta, error%beta2 /)
        endif

        if (nparam == 16) then
            values_vec = (/ error%bgauss, error%thetabd, error%chibd, error%bgauss2, error%thetabd2, error%chibd2,&
                error%height, error%dtau, error%dtau2, error%vdopp, error%vdopp2, error%damping, &
                error%vmacro, error%vmacro2, error%beta, error%beta2 /)
        endif

        if (nparam == 17) then
            values_vec = (/ error%bgauss, error%thetabd, error%chibd, error%bgauss2, error%thetabd2, error%chibd2,&
                error%height, error%dtau, error%dtau2, error%vdopp, error%vdopp2, error%damping, &
                error%vmacro, error%vmacro2, error%beta, error%ff, error%beta2 /)
        endif
        
        call check( nf90_put_var(in_fixed%error_id, in_fixed%map_error_id, values_vec, start=start, count=count) )

        status = nf90_sync(in_fixed%syn_id)
        status = nf90_sync(in_fixed%par_id)
        status = nf90_sync(in_fixed%error_id)
        
        deallocate(values_vec, start, count)
#endif
        
    end subroutine write_results

!------------------------------------------------------------
! Read the file with the parameters to invert
!------------------------------------------------------------
    subroutine read_parameters_to_invert(file,in_params,in_inversion,in_fixed)
    character(len=120) :: file
    type(variable_parameters) :: in_params
    type(type_inversion) :: in_inversion
    type(fixed_parameters) :: in_fixed
    integer :: i, j
        
! Number of parameters of the inversion     
        in_params%n_total = 18

        allocate(in_params%inverted(in_params%n_total))
        in_params%inverted = 0
                
        open(unit=12,file=file,action='read',status='old')
        call lb(12,3)
        call lb(12,2)
        read(12,*) in_inversion%iter_max
        call lb(12,2)
        read(12,*) in_inversion%n_cycles

        allocate(in_inversion%cycles(in_params%n_total,in_inversion%n_cycles))
        allocate(in_inversion%stokes_weights(0:3,in_inversion%n_cycles))

        do i = 1, in_params%n_total
            call lb(12,2)
            read(12,*) (in_inversion%cycles(i,j),j=1,in_inversion%n_cycles)

! Test whether the damping is put in automatic and inverted. If this is the case, make it not invertible
            if (i == 8) then
                if (in_fixed%damping_treatment == 1) then
                    write(*,*) 'Damping is automatic. Unsetting inversion switch'
                    in_inversion%cycles(i,:) = 0
                endif
            endif
        enddo

        do i = 0, 3
            call lb(12,2)
            read(12,*) (in_inversion%stokes_weights(i,j),j=1,in_inversion%n_cycles)
        enddo
        
        call lb(12,2)
        allocate(in_inversion%algorithm(in_inversion%n_cycles))
        read(12,*) (in_inversion%algorithm(j),j=1,in_inversion%n_cycles)

                
        close(12)

        in_params%inverted = in_inversion%cycles(:,1)
    
        in_params%n_inverted = sum(params%inverted)

    end subroutine read_parameters_to_invert
    
!------------------------------------------------------------
! Write a file with the final parameters so that it can be used for initializing again the inversion
!------------------------------------------------------------
    subroutine write_experiment(in_params,in_fixed,err_params)
    character(len=120) :: tmp
    type(variable_parameters) :: in_params, err_params
    type(fixed_parameters) :: in_fixed
    integer :: i
    
        open(unit=13,file=output_final_parameters,action='write',status='replace')

        write(13,FMT='(A)') '***********************************************'
        write(13,FMT='(A)') '* File defining the specific experiment to solve'
        write(13,FMT='(A)') '***********************************************'
        write(13,*)
        write(13,FMT='(A)') '# Include stimulated emission (0-> no, 1-> yes)'
        write(13,FMT='(I1)') isti
        write(13,*)
        
        write(13,FMT='(A)') '# Include magnetic field (0-> no, 1-> yes)'
        write(13,FMT='(I1)') imag
        write(13,*)

        write(13,FMT='(A)') '# Include depolarization rates (0-> no, 1-> yes)'
        write(13,FMT='(I1)') idep
        write(13,*)

        write(13,FMT='(A)') '# Value of delta if depolarization rates are included (not used if the previous value is 0)'
        write(13,FMT='(E10.3)') in_params%delta_collision
        write(13,*)
        
        write(13,FMT='(A)') '# Include Paschen-Back effect (0-> no, 1-> yes)'
        write(13,FMT='(I1)') use_paschen_back
        write(13,*)

        write(13,FMT='(A)') '# Number of slabs (1-> 1 slab, 2-> 2 slabs with same B, 3-> 2 slabs with different B (add field below))'
        write(13,FMT='(I1)') in_params%nslabs
        write(13,*)
        
        if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
            write(13,FMT='(A)') '# Magnetic field strength [G], thetaB [degrees], chiB [degrees]'
            write(13,FMT='(F9.4,3X,F7.2,3X,F7.2)') in_params%bgauss, in_params%thetabd, in_params%chibd
            write(13,*)
        endif
        if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(A)') '# Magnetic field strength [G], thetaB [degrees], chiB [degrees], B2 [G], thetaB2 [degrees], chiB2 [degrees]'
            write(13,FMT='(F9.4,3X,F7.2,3X,F7.2,3X,F9.4,3X,F7.2,3X,F7.2)') in_params%bgauss, in_params%thetabd, in_params%chibd, in_params%bgauss2, in_params%thetabd2, in_params%chibd2
            write(13,*)
        endif
        
        write(13,FMT='(A)') '# Apparent height of the He I atoms in arcsec'
        write(13,FMT='(F7.3)') in_params%height
        write(13,*)
        
        write(13,FMT='(A)') '# Optical depth of the slab in the maximum of I (slab) or strength of the line (ME)'
        if (in_params%nslabs == 1) then
            write(13,FMT='(F8.4)') in_params%dtau
        endif
        if (in_params%nslabs == 2 .or. in_params%nslabs == 3) then
            write(13,FMT='(F8.4,2X,F8.4)') in_params%dtau, in_params%dtau2
        endif
        if (in_params%nslabs == -2) then
            write(13,FMT='(F8.4,2X,F8.4)') in_params%dtau, in_params%dtau2, in_params%ff
        endif
        write(13,*)
                
        write(13,FMT='(A)') '# Source function gradient (only ME)'
        if (in_params%nslabs == 2 .or. in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(F8.4)') in_params%beta, in_params%beta2
        else
            write(13,FMT='(F8.4)') in_params%beta
        endif
        write(13,*)
        
        write(13,FMT='(A)') '# Boundary Stokes parameters (I0,Q0,U0,V0)'
        write(13,FMT='(4(E10.3,2X))') (in_fixed%Stokes_incident(i),i=0,3)
        write(13,*)
        
        write(13,FMT='(A)') '# Transition where to compute the emission'
        write(13,FMT='(I1)') in_fixed%nemiss
        write(13,*)
        
        write(13,FMT='(A)') '# Use atomic polarization? (0-> no, 1-> yes)'
        write(13,FMT='(I1)') in_fixed%use_atomic_pol
        write(13,*)
        
        write(13,FMT='(A)') '# Observation angle with respect to the vertical theta,chi,gamma [degrees]'
        write(13,FMT='(3(F7.2,2X))') in_fixed%thetad, in_fixed%chid, in_fixed%gammad
        write(13,*)
        
        write(13,FMT='(A)') '# Wavelength axis: minimum, maximum and number of grid points'
        write(13,*) in_fixed%omin, in_fixed%omax, in_fixed%no
        write(13,FMT='(3(F7.2,2X))')

        if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
            write(13,FMT='(A)') '# Line wavelength [A], Doppler velocity [km/s] and damping [a]'
            write(13,FMT='(F10.4,3X,F6.3,3X,F7.4)') in_fixed%wl, in_params%vdopp, in_params%damping
            write(13,*)
        endif
        if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(A)') '# Line wavelength [A], Doppler velocity [km/s], Doppler velocity 2 [km/s] and damping [a]'
            write(13,FMT='(F10.4,3X,F6.3,3X,F6.3,3X,F7.4)') in_fixed%wl, in_params%vdopp, in_params%vdopp2, in_params%damping
            write(13,*)
        endif
        
        write(13,FMT='(A)') '# Macroscopic velocity [km/s] (>0 is a redshift)'
        if (in_params%nslabs == 1) then
            write(13,FMT='(F7.3)') in_params%vmacro
        else
            write(13,FMT='(F7.3,2X,F7.3)') in_params%vmacro, in_params%vmacro2
        endif
        write(13,*)
        
        write(13,FMT='(A)') '# Include magneto-optical effects in the RT'
        write(13,FMT='(I1)') use_mag_opt_RT
        write(13,*)
        
        write(13,FMT='(A)') '# Include stimulated emission in the RT'
        write(13,FMT='(I1)') use_stim_emission_RT
        write(13,*)
        
        close(13)



        open(unit=13,file=output_error_parameters,action='write',status='replace')

        write(13,FMT='(A)') '***********************************************'
        write(13,FMT='(A)') '* File defining the specific experiment to solve'
        write(13,FMT='(A)') '***********************************************'
        write(13,*)
        write(13,FMT='(A)') '# Include stimulated emission (0-> no, 1-> yes)'
        write(13,FMT='(I1)') isti
        write(13,*)
        
        write(13,FMT='(A)') '# Include magnetic field (0-> no, 1-> yes)'
        write(13,FMT='(I1)') imag
        write(13,*)

        write(13,FMT='(A)') '# Include depolarization rates (0-> no, 1-> yes)'
        write(13,FMT='(I1)') idep
        write(13,*)

        write(13,FMT='(A)') '# Value of delta if depolarization rates are included (not used if the previous value is 0)'
        write(13,FMT='(E10.3)') err_params%delta_collision
        write(13,*)
        
        write(13,FMT='(A)') '# Include Paschen-Back effect (0-> no, 1-> yes)'
        write(13,FMT='(I1)') use_paschen_back
        write(13,*)

        write(13,FMT='(A)') '# Number of slabs (1-> 1 slab, 2-> 2 slabs with same B, 3-> 2 slabs with different B (add field below))'
        write(13,FMT='(I1)') err_params%nslabs
        write(13,*)
        
        if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
            write(13,FMT='(A)') '# Magnetic field strength [G], thetaB [degrees], chiB [degrees]'
            write(13,FMT='(F9.4,3X,F7.2,3X,F7.2)') err_params%bgauss, err_params%thetabd, err_params%chibd
            write(13,*)
        endif
        if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(A)') '# Magnetic field strength [G], thetaB [degrees], chiB [degrees], B2 [G], thetaB2 [degrees], chiB2 [degrees]'
            write(13,FMT='(F9.4,3X,F7.2,3X,F7.2,3X,F9.4,3X,F7.2,3X,F7.2)') err_params%bgauss, err_params%thetabd, err_params%chibd, err_params%bgauss2, err_params%thetabd2, err_params%chibd2
            write(13,*)
        endif
        
        write(13,FMT='(A)') '# Apparent height of the He I atoms in arcsec'
        write(13,FMT='(F7.3)') err_params%height
        write(13,*)
        
        write(13,FMT='(A)') '# Optical depth of the slab in the maximum of I (slab) or strength of the line (ME)'
        if (in_params%nslabs == 1) then
            write(13,FMT='(F8.4)') err_params%dtau
        endif
        if (in_params%nslabs == 2 .or. in_params%nslabs == 3) then
            write(13,FMT='(F8.4,2X,F8.4)') err_params%dtau, err_params%dtau2
        endif
        if (in_params%nslabs == -2) then
            write(13,FMT='(F8.4,2X,F8.4)') err_params%dtau, err_params%dtau2, err_params%ff
        endif
        write(13,*)
                
        write(13,FMT='(A)') '# Source function gradient (only ME)'
        if (in_params%nslabs == 2 .or. in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(F8.4)') err_params%beta, err_params%beta2
        else
            write(13,FMT='(F8.4)') err_params%beta
        endif
        write(13,*)
        
        write(13,FMT='(A)') '# Boundary Stokes parameters (I0,Q0,U0,V0)'
        write(13,FMT='(4(E10.3,2X))') (in_fixed%Stokes_incident(i),i=0,3)
        write(13,*)
        
        write(13,FMT='(A)') '# Transition where to compute the emission'
        write(13,FMT='(I1)') in_fixed%nemiss
        write(13,*)
        
        write(13,FMT='(A)') '# Use atomic polarization? (0-> no, 1-> yes)'
        write(13,FMT='(I1)') in_fixed%use_atomic_pol
        write(13,*)
        
        write(13,FMT='(A)') '# Observation angle with respect to the vertical theta,chi,gamma [degrees]'
        write(13,FMT='(3(F7.2,2X))') in_fixed%thetad, in_fixed%chid, in_fixed%gammad
        write(13,*)
        
        write(13,FMT='(A)') '# Wavelength axis: minimum, maximum and number of grid points'
        write(13,*) in_fixed%omin, in_fixed%omax, in_fixed%no
        write(13,FMT='(3(F7.2,2X))')

        if (in_params%nslabs == 1 .or. in_params%nslabs == 2) then
            write(13,FMT='(A)') '# Line wavelength [A], Doppler velocity [km/s] and damping [a]'
            write(13,FMT='(F10.4,3X,F6.3,3X,F7.4)') in_fixed%wl, err_params%vdopp, err_params%damping
            write(13,*)
        endif
        if (in_params%nslabs == 3 .or. in_params%nslabs == -2) then
            write(13,FMT='(A)') '# Line wavelength [A], Doppler velocity [km/s], Doppler velocity 2 [km/s] and damping [a]'
            write(13,FMT='(F10.4,3X,F6.3,3X,F6.3,3X,F7.4)') in_fixed%wl, err_params%vdopp, err_params%vdopp2, err_params%damping
            write(13,*)
        endif
        
        write(13,FMT='(A)') '# Macroscopic velocity [km/s] (>0 is a redshift)'
        if (in_params%nslabs == 1) then
            write(13,FMT='(F7.3)') err_params%vmacro
        else
            write(13,FMT='(F7.3,2X,F7.3)') err_params%vmacro, err_params%vmacro2
        endif
        write(13,*)
        
        write(13,FMT='(A)') '# Include magneto-optical effects in the RT'
        write(13,FMT='(I1)') use_mag_opt_RT
        write(13,*)
        
        write(13,FMT='(A)') '# Include stimulated emission in the RT'
        write(13,FMT='(I1)') use_stim_emission_RT
        write(13,*)
        
        close(13)
    
    end subroutine write_experiment 
    
!------------------------------------------------------------
! Save the inverted profiles
!------------------------------------------------------------
    subroutine write_final_profiles(file,in_observation,in_inversion)
    character(len=120) :: file
    type(type_observation) :: in_observation
    type(type_inversion) :: in_inversion
    integer :: i, j
    
        open(unit=12,file=file,action='write',status='replace')
        
        write(12,*) in_observation%n, in_inversion%chisq
        
        do i = 1, in_observation%n
            write(12,FMT='(5(1X,E15.8))') in_observation%wl(i), &
                (in_inversion%stokes_unperturbed(j,i),j=0,3)
        enddo
        
        close(12)
    end subroutine write_final_profiles 

!------------------------------------------------------------
! Print the parameters depending on the synthesis mode
!------------------------------------------------------------
    subroutine print_parameters(params,text,header)
    type(variable_parameters) :: params
    character(len=20) :: text
    logical :: header

        if (params%nslabs == 1) then
! 1 component
! EMISSION
            if (synthesis_mode == 0) then
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       D^(2)    vmacro      a        h'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       vmacro       a       h'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%vmacro, params%damping, params%height
                endif               
            endif
    ! CONSTANT SLAB 
            if (synthesis_mode == 1 .or. synthesis_mode == 3 .or. synthesis_mode == 4 .or. &
                synthesis_mode == 5) then
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         tau       D^(2)     vmacro       a       h'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,3X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%dtau, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         beta      tau       vmacro       a        h'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%beta, params%dtau, params%vmacro, params%damping, params%height
                endif               
            endif
        endif

        if (params%nslabs == 2) then
! 2 components
! EMISSION
            if (synthesis_mode == 0) then                   
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       D^(2)    vmacro      a        h      v_th2     beta      beta2'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height, params%vmacro2, params%beta, params%beta2
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th       vmacro       a       h     v_th2     beta      beta2'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%vmacro, params%damping, params%height, params%vmacro2, params%beta, params%beta2
                endif               
            endif
    ! CONSTANT SLAB 
            if (synthesis_mode == 1 .or. synthesis_mode == 3 .or. synthesis_mode == 4 .or. &
                synthesis_mode == 5) then
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         tau       D^(2)     vmacro       a       h        tau2     vmacro2     beta      beta2'
                    write(*,FMT='(A,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,3X,E9.2,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4)') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%dtau, 10.d0**params%delta_collision, params%vmacro, params%damping, params%height, params%dtau2, params%vmacro2, params%beta, params%beta2                   
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        v_th         tau       vmacro       a        h        tau2     vmacro2     beta      beta2'
                    write(*,FMT='(A,12(F9.4,2X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%vdopp, params%dtau, params%vmacro, params%damping, params%height, params%dtau2, params%vmacro2, params%beta, params%beta2
                    
                endif               
            endif
        endif

        if (params%nslabs == 3) then
! 2 components with different fields
! EMISSION
            if (synthesis_mode == 0) then                   
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB      B2      thetaB2     chiB2      v_th      v_th2      D^(2)    vmacro     vmacro2    a        h     beta      beta2'
                    write(*,FMT='(A,6(F9.4,1X),E9.2,1X,5(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, 10.d0**params%delta_collision, &
                        params%vmacro, params%vmacro2, params%damping, params%height, params%beta, params%beta2
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB      B2      thetaB2     chiB2      v_th     v_th2     vmacro      vmacro2     a       h     beta      beta2'
                    write(*,FMT='(A,14(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%vmacro, params%vmacro2, params%damping, params%height, params%beta, params%beta2
                endif               
            endif
    ! CONSTANT SLAB 
            if (synthesis_mode == 1 .or. synthesis_mode == 3 .or. synthesis_mode == 4 .or. synthesis_mode == 5) then
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB       B2      thetaB2     chiB2      v_th      v_th2     tau     tau2     D^(2)     vmacro      vmacro2     a       h     beta      beta2'
                    write(*,FMT='(A,10(F9.4,1X),E9.2,1X,4(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%dtau, params%dtau2, &
                        10.d0**params%delta_collision, params%vmacro, params%vmacro2, params%damping, params%height, params%beta, params%beta2
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        B2     thetaB2    chiB2     v_th      v_th2     tau      tau2      vmacro    vmacro2       a        h     beta      beta2'
                    write(*,FMT='(A,16(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%dtau, params%dtau2, params%vmacro, params%vmacro2,&
                        params%damping, params%height, params%beta, params%beta2
                endif               
            endif
        endif
        
        if (params%nslabs == -2) then
! 2 components with different fields
! EMISSION
            if (synthesis_mode == 0) then                   
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB      B2      thetaB2     chiB2      v_th      v_th2      D^(2)    vmacro     vmacro2   ff    a        h     beta      beta2'
                    write(*,FMT='(A,6(F9.4,1X),E9.2,1X,7(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, 10.d0**params%delta_collision, &
                        params%vmacro, params%vmacro2, params%ff, params%damping, params%height, params%beta, params%beta2
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB      B2      thetaB2     chiB2      v_th     v_th2     vmacro      vmacro2     ff    a       h     beta      beta2'
                    write(*,FMT='(A,15(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%vmacro, params%vmacro2, params%ff, params%damping, params%height, params%beta, params%beta2
                endif               
            endif
! CONSTANT SLAB 
            if (synthesis_mode == 1 .or. synthesis_mode == 3 .or. synthesis_mode == 4 .or. &
                synthesis_mode == 5) then
                if (params%inverted(6) == 1) then
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB       B2      thetaB2     chiB2      v_th      v_th2     tau     tau2     D^(2)     vmacro      vmacro2      ff      a        h     beta      beta2'
                    write(*,FMT='(A,10(F9.4,1X),E9.2,1X,7(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%dtau, params%dtau2, &
                        10.d0**params%delta_collision, params%vmacro, params%vmacro2, params%ff, params%damping, params%height, params%beta, params%beta2
                else
                    if (header) write(*,FMT='(A)') '                        B        thetaB     chiB        B2     thetaB2    chiB2     v_th      v_th2     tau      tau2      vmacro    vmacro2      ff       a         h     beta      beta2'
                    write(*,FMT='(A,17(F9.4,1X))') text, params%bgauss, params%thetabd, &
                        params%chibd, params%bgauss2, params%thetabd2, params%chibd2, params%vdopp, params%vdopp2, params%dtau, params%dtau2, params%vmacro, params%vmacro2,&
                        params%ff, params%damping, params%height, params%beta, params%beta2
                endif               
            endif
        endif
    
        
    end subroutine print_parameters
end module io
