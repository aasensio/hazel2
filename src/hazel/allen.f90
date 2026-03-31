module allen
use vars
use maths
implicit none

contains
!------------------------------------------------------------
! Read Allen's data
!------------------------------------------------------------
	subroutine read_allen_data
	integer :: i, j
			
		! open(unit=12, file='CLV/ic.dat', action='read', status='old')
		! do i = 1, 43
		! 	read(12, *) (allen_ic(i, j), j=1, 2)
		! enddo
		! close(12)

		if (allocated(allen_ic)) deallocate(allen_ic)
		if (allocated(allen_cl)) deallocate(allen_cl)
		allocate(allen_ic(43,2))
		allocate(allen_cl(22,3))
		
		allen_cl(:,1) = (/0.200000,0.220000,0.245000,0.265000,0.280000,0.300000,0.320000,0.350000,0.370000,0.380000,&
			0.400000,0.450000,0.500000,0.550000,0.600000,0.800000,1.00000,1.50000,2.00000,3.00000,5.00000,10.0000/)
		allen_cl(:,2) = (/0.120000,-1.30000,-0.100000,-0.100000,0.380000,0.740000,0.880000,0.980000,1.03000,0.920000,&
			0.910000,0.990000,0.970000,0.930000,0.880000,0.730000,0.640000,0.570000,0.480000,0.350000,0.220000,0.150000/)
		allen_cl(:,3) = (/0.330000,1.60000,0.850000,0.900000,0.570000,0.200000,0.0300000,-0.100000,-0.160000,-0.0500000,&
			-0.0500000,-0.170000,-0.220000,-0.230000,-0.230000,-0.220000,-0.200000,-0.210000,-0.180000,-0.120000,-0.0700000,-0.0700000/)
			
		allen_ic(:,1) = (/0.200000,0.220000,0.240000,0.260000,0.280000,0.300000,0.320000,0.340000,0.360000,0.370000,0.380000,0.390000,&
			0.400000,0.410000,0.420000,0.430000,0.440000,0.450000,0.460000,0.480000,0.500000,0.550000,0.600000,0.650000,0.700000,0.750000,&
			0.800000,0.900000,1.00000,1.10000,1.20000,1.40000,1.60000,1.80000,2.00000,2.50000,3.00000,4.00000,5.00000,6.00000,8.00000,10.0000,12.0000/)
			
		allen_ic(:,2) = (/0.0600000,0.210000,0.290000,0.600000,1.30000,2.45000,3.25000,3.77000,4.13000,4.23000,4.63000,4.95000,5.15000,5.26000,5.28000,&
			5.24000,5.19000,5.10000,5.00000,4.79000,4.55000,4.02000,3.52000,3.06000,2.69000,2.28000,2.03000,1.57000,1.26000,1.01000,0.810000,0.530000,0.360000,&
			0.238000,0.160000,0.0780000,0.0410000,0.0142000,0.00620000,0.00320000,0.000950000,0.000350000,0.000180000/)

! Units conversion
		allen_ic(:, 2) = allen_ic(:, 2) * allen_ic(:, 1)**2 / (PC*1.d4) ! i_lambda to i_nu
		allen_ic(:, 1) = 1d4 * allen_ic(:, 1)  ! wavelength in angstrom

! Reading data file with coefficients for centre-to-limb variation
!      first column: wavelength in microns
!      2nd & 3rd columns: u & v coefficients where
!         i(mu)/i(0) = 1 - u - v + u \mu + v \mu^2
		! open(unit=12, file='CLV/cl.dat', action='read', status='old')
  !     do i = 1, 22
		! 	read(12, *) (allen_cl(i, j), j=1, 3)
  !     enddo
  !     close(12)

      allen_cl(:, 1) = 1d4 * allen_cl(:, 1)  ! wavelength in angstrom
		
	end subroutine read_allen_data
	
!------------------------------------------------------------
! Return the nbar parameter from Allen's data
!------------------------------------------------------------
	function nbar_allen(lmb, in_fixed, in_params, reduction_factor)
	type(fixed_parameters) :: in_fixed
	type(variable_parameters) :: in_params
	real(kind=8) :: lmb, height, nbar_allen, l(1), t(1), delta, reduction_factor
	integer :: i
	real(kind=8) :: u1, u2, I0, sg, cg, a0, a1, a2, b0, b1, b2, J, K

		l(1) = lmb
		call lin_interpol(allen_ic(:,1), allen_ic(:,2), l, t)
		I0 = t(1)
		call lin_interpol(allen_cl(:,1), allen_cl(:,2), l, t)
		u1 = t(1)
		call lin_interpol(allen_cl(:,1), allen_cl(:,3), l, t)
		u2 = t(1)
		
! Take into account the observation angle (scattering angle) to calculate the height from the apparent height
! Only do this if we use the apparent height
		if (in_params%height < 0.d0) then
			delta = abs(90.d0 - in_fixed%thetad)
			height = (RSUN + in_params%height) / (cos(delta*PI/180.d0)) - RSUN
		else
			height = abs(in_params%height)
		endif

		if (height /= 0.d0) then
      	sg = Rsun / (height+Rsun)
      	cg = sqrt(1.d0-sg**2)

      	a0 = 1.d0 - cg
      	a1 = cg - 0.5d0 - 0.5d0*cg**2/sg*log((1.d0+sg)/cg)
      	a2 = (cg+2.d0)*(cg-1.d0) / (3.d0*(cg+1.d0))
		else
			a0 = 1.d0
			a1 = -0.5d0
			a2 = -2.d0 / 3.d0
		endif
        
		J = 0.5d0 * I0 * (a0 + a1*u1 + a2*u2)

		nbar_allen = 1d10*(lmb**3*1d-24/(2d0*ph*pc)) * J * reduction_factor
		
	end function nbar_allen
	
!------------------------------------------------------------
! Return the anisotropy parameter from Allen's data
! Reduction factor affects only the calculation of J00
!------------------------------------------------------------
	function omega_allen(lmb, in_fixed, in_params, reduction_factor)
	type(fixed_parameters) :: in_fixed
	type(variable_parameters) :: in_params
	real(kind=8) :: lmb, height, omega_allen, l(1), t(1), delta, reduction_factor
	integer :: i
	real(kind=8) :: u1, u2, I0, sg, cg, a0, a1, a2, b0, b1, b2, J, K

		l(1) = lmb
		call lin_interpol(allen_ic(:,1), allen_ic(:,2), l, t)
		I0 = t(1)
		call lin_interpol(allen_cl(:,1), allen_cl(:,2), l, t)
		u1 = t(1)
		call lin_interpol(allen_cl(:,1), allen_cl(:,3), l, t)
		u2 = t(1)
		
	
! Take into account the observation angle (scattering angle) to calculate the height from the apparent height
! Only do this if we use the apparent height
		if (in_params%height < 0.d0) then
			delta = abs(90.d0 - in_fixed%thetad)
			height = (RSUN + in_params%height) / (cos(delta*PI/180.d0)) - RSUN
		else
			height = abs(in_params%height)
		endif
		
      if (height /= 0) then
			sg = Rsun / (height+Rsun)
      	cg = sqrt(1.d0-sg**2)

      	a0 = 1.d0 - cg
      	a1 = cg - 0.5d0 - 0.5d0*cg**2/sg*log((1.d0+sg)/cg)
      	a2 = (cg+2.d0)*(cg-1.d0) / (3.d0*(cg+1.d0))

      	b0 = (1.d0-cg**3) / 3.d0
      	b1 = (8.d0*cg**3-3.d0*cg**2-2.d0) / 24.d0 - cg**4 / (8.d0*sg) * log((1.d0+sg)/cg)
      	b2 = (cg-1.d0)*(3.d0*cg**3+6.d0*cg**2+4.d0*cg+2.d0) / (15.d0*(cg+1.d0))
		else
			a0 = 1.d0
			a1 = -0.5d0
			a2 = -2.d0 / 3.d0
			
			b0 = 1.d0 / 3.d0
			b1 = -1.d0 / 12.d0
			b2 = -2.d0 / 15.d0
		endif
        
		J = 0.5d0 * I0 * (a0 + a1*u1 + a2*u2)
		K = 0.5d0 * I0 * (b0 + b1*u1 + b2*u2)

		omega_allen = (3.d0*K-J) / (2.d0*J) * reduction_factor
		
	end function omega_allen	

!------------------------------------------------------------
! Return the continuum intensity from Allen's data
!------------------------------------------------------------
	function I0_allen(lmb, mu)
	real(kind=8) :: I0_allen
	real(kind=8) :: lmb, mu, l(1), t(1), delta, reduction_factor
	integer :: i
	real(kind=8) :: u1, u2, I0, sg, cg, a0, a1, a2, b0, b1, b2, J, K

		l(1) = lmb
		call lin_interpol(allen_ic(:,1), allen_ic(:,2), l, t)
		I0 = t(1) * 1.d10
		call lin_interpol(allen_cl(:,1), allen_cl(:,2), l, t)
		u1 = t(1)
		call lin_interpol(allen_cl(:,1), allen_cl(:,3), l, t)
		u2 = t(1)
		
		I0_allen = I0*(1.d0 - u1 - u2 + u1 * mu + u2 * mu**2)
		
	end function I0_allen
	
end module allen
