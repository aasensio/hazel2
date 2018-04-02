c rutina delta1

	subroutine delta1(ntau,tk,ts,us,oo,dab,rt)

	include 'PARAMETER'  !por kt
	implicit real*4 (a-h,o-z)
	parameter (kt16=16*kt)

	real*4 tk(kt16),ts(kt),us(4,kt),oo(4,4,kt)
        real*4 dab(kt16),rt(4,kt)
	real*4 trhs(4,kt)

c calculamos trhs=(-[delta k]* (i-s) + [k]*(delta s) )
	do k=1,ntau
	   do i=1,4
	      trhs(i,k)=0.
	      do j=1,4
	         k1=i+4*j+16*k-20	              !(i,j,k)
	         trhs(i,k)=trhs(i,k)-tk(k1)*us(j,k)   !-[deltak]*(i-s)
	      end do
	      k2=i+16*k-16			      !(i,1,k)
	      trhs(i,k)= trhs(i,k) + dab(k2)*ts(k)
	   end do
	end do

c calculamos la f. respuesta rt(k)=[o](k)*j(k)
	do k=1,ntau-1
	   do i=1,4
	      a=0.
	      do j=1,4
	         a=a+oo(i,j,k)*trhs(j,k)
	      end do
	      rt(i,k)=-a    !el signo viene de que integro al reves
	   end do
	end do
  	do i=1,4
  	   rt(i,ntau)=0.0		!rt(i,ntau-1)
  	end do

	return
	end
