      subroutine matinx ( a )

      
c	'exact' inversion of 4 X 4 matrix

      implicit real*4 ( a-h, o-z )

      dimension a ( 4 , 4 ) , b ( 4 , 4 )
c      dimension c(4,4)
      integer i , j

      integer :: error_code
	common/Error/error_code

      absmax = 0.

      do 5010 i = 1 , 4
      do 5010 j = 1 , 4
c          c(i,j)=a(i,j) 
	  if ( abs ( a ( i , j ) ) .gt. absmax ) absmax =abs( a ( i , j ))
5010  continue

c      if ( absmax .eq. 0. ) stop 'singularity problem in matinx'

      if ( absmax .eq. 0. ) then
            error_code = 1
            a = 1.0
            return
      endif


      fabsmx = 1.d00 / absmax


      do 5020 i = 1 , 4
      do 5020 j = 1 , 4
	  a ( i , j ) = a ( i , j ) * fabsmx
5020  continue



      b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)
     1 + a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)
     2 - a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
      b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)
     1 + a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)
     2 - a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
      b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)
     1 + a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)
     2 - a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
      b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)
     1 + a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)
     2 - a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
      b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)
     1 + a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)
     2 - a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
      b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)
     1 + a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)
     2 - a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
      b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)
     1 + a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)
     2 - a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
      b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)
     1 + a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)
     2 - a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
      b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)
     1 + a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)
     2 - a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
      b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)
     1 + a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)
     2 - a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
      b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)
     1 + a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)
     2 - a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
      b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)
     1 + a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)
     2 - a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
      b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)
     1 + a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)
     2 - a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
      b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)
     1 + a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)
     2 - a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
      b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)
     1 + a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)
     2 - a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
      b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)
     1 + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)
     2 - a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

      det = a ( 1 , 1 ) * b ( 1 , 1 ) + a ( 1 , 2 ) * b ( 2 , 1 )
     1 + a ( 1 , 3 ) * b ( 3 , 1 ) + a ( 1 , 4 ) * b ( 4 , 1 )
c	if(abs(det).lt.1.e-2)then
c		print*,det
c                do i=1,4
c                   do j=1,4
c                       print*,'c(',i,',',j,')=',sngl(c(i,j))  
c                   end do
c                end do
c	end if

c      print*,'en matinx (4)',det

      fdeta = fabsmx / det

c      print*,'en matinx (5)',fdeta


      do 5000 i = 1 , 4
      do 5000 j = 1 , 4
          a ( i , j ) = b ( i , j ) * fdeta
5000  continue
c      	if(abs(det).lt.1.e-2)then
c          do i=1,4
c             do j=1,4
c                sum=0.
c                do l=1,4
c                    sum=sum+a(i,l)*c(l,j)
c                end do
c                print*,'c*inversa(',i,',',j,')=',sum
c             end do
c           end do
c         end if 

c      print*,'en matinx (6)',fdeta
   
      return
      end

c______________________________________________________________________________
      subroutine matinx4 ( b )
c______________________________________________________________________________

      implicit real*4 ( a-h, o-z )

      dimension b ( 4 , 4 ) 
	
	q=b(1,2)
	u=b(1,3)
	v=b(1,4)
	r=b(2,3)
	s=b(4,2)
	t=b(3,4)
	q2=q*q
	u2=u*u
	v2=v*v
	r2=r*r
	s2=s*s
	t2=t*t
	a=q*t+r*v+s*u
	a2=a*a
	taq=t*a+q
	qat=q*a-t
	sau=s*a+u
	uas=u*a-s
	rav=r*a+v
	var=v*a-r
	ur=u*r-s*v
	uv=u*v+r*s
	qr=q*r-t*v
	qs=q*s-t*u
	qu=q*u+t*s
	qv=q*v+t*r

	det=1.+t2+r2+s2-q2-u2-v2-a2
	
	b(1,1)=(1.+t2+r2+s2)/det
	b(2,1)=(ur-taq)/det
	b(3,1)=(-qr-sau)/det
	b(4,1)=(qs-rav)/det
	b(1,2)=(-ur-taq)/det
	b(2,2)=(1.+t2-u2-v2)/det
	b(3,2)=(qu-var)/det
        b(4,2)=(qv+uas)/det
	b(1,3)=(qr-sau)/det
	b(2,3)=(qu+var)/det
	b(3,3)=(1.+s2-q2-v2)/det
	b(4,3)=(uv-qat)/det
	b(1,4)=(-qs-rav)/det
        b(2,4)=(qv-uas)/det
	b(3,4)=(uv+qat)/det
	b(4,4)=(1.+r2-q2-u2)/det

      return
      end
