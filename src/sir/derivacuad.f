	subroutine derivacuad(s,ds,n,x)

	integer n
	real*4 s(*), ds(*),x(*),x1,x2,x3,d1,d2,d3

	ds(1)=(s(2)-s(1))/(x(2)-x(1))
        do i=2,n-1
           x1=x(i-1)
           x2=x(i)
           x3=x(i+1)
           d1=1./(x2-x3)
           d2=1./(x1-x2)
           d3=1./(x1-x3)
           ds(i)=s(i+1)*(d3-d1)+s(i)*(d1-d2)+s(i-1)*(d2-d3)
        end do
	ds(n)=(s(n)-s(n-1))/(x(n)-x(n-1))


	return
	end
