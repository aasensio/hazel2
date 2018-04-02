	subroutine deriva(s,ds,n,nn,TOTAL)

	integer n
	real*4 s(*), ds(*)
	real*4 TOTAL(nn,nn)


c	do i=1,n
c	   ds(i)=0.
c	   do j=1,n
c	      ds(i) = ds(i) + TOTAL(i,j)*s(j)
c	   enddo
c	enddo

	do i=1,n
	   dss=0.
	   do j=1,n
	      dss = dss + TOTAL(i,j)*s(j)
	   enddo
           ds(i)=dss
	enddo

	return
	end
