c rnorma rutina que reescala las funciones con el paso de la integral
c
	subroutine rnorma(ntau,cont,r)

	include 'PARAMETER'   !para kt
	parameter (aln10=2.3025851)
	implicit real*4 (a-h,o-z)
	real*4 tau(kt),taue(kt),r(*)
	real*4 deltae(kt),deltai(kt),delt2i(kt)

        common/segunda/tau,taue,deltae,deltai,delt2i


	paso= (tau(2) - tau(1) ) *aln10/cont
	do i=1,ntau
	   r(i)=r(i)*taue(i)*paso
	end do

	return
	end
