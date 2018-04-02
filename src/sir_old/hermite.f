	subroutine hermite(s,ds,etall,n,svec,nn,bt,tk,pk,hk,
     &           vk,gk,fk,mk,rt,rp,rh,rv,rg,rf,rm,mnodos)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER' !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

	integer n,nn
	real*4 svec(4)
        integer*4 mnodos(*)

c internas
        real*4 sol(4,kt),tau(kt),taue(kt),s0(kt)
        real*4 deltae(kt),deltai(kt),delt2i(kt)
	real*4 ds(kt), detal(4,4,kt), etal(4,4,kt), etall(4,4,kt)
	real*4 etal2(4,4,kt)
	real*4 s(kt), den(4), num(4,4), den2(4,4)
	real*4 OO(4,4,kt)
        real*4 tk(kt16),pk(kt16),hk(kt16),vk(kt16),gk(kt16),fk(kt16),mk(kt16)
	real*4 bt(kt)
	real*4 rt(4,kt),rp(4,kt),rh(4,kt),rv(4,kt),rg(4,kt),rf(4,kt),rm(4,kt)

        common/segunda/tau,taue,deltae,deltai,delt2i


c se calcula donde poner el contorno

        icontorno=1
	do i=1,n-1

	   ti=abs(etall(1,1,i))
           tq=abs(etall(1,2,i))
           tu=abs(etall(1,3,i))
           tv=abs(etall(1,4,i))
c	   print*,'ti,tq,tu,tv=',i,ti,tq,tu,tv


           taunu=max(ti,tq,tu,tv)*deltae(i)
           if(taunu.gt.5.)icontorno=i

        end do
c	icontorno=1


c se calcula la derivada de la matriz de absorcion

	do i=1,n
	   xx = taue(i)*dl10
c	   print*,'tau,taue,deltae,deltai,delt2i=',tau(i),taue(i),deltae(i),deltai(i),delt2i(i)
	   do jj=1,4
	      do kk=1,4
	         etal(jj,kk,i) = etall(jj,kk,i)*xx
c		 print*,'etal,jj,kk,i=',jj,kk,i,etal(jj,kk,i),xx
	      enddo
	   enddo
  	end do 	
	call deriva4cuad(etal,detal,n,tau)


c se integra
        call inicializa(icontorno,etal,s,ds,oo,sol,ierror)
	if(ierror.eq.1)return
c	ic=icontorno
c	do while(ierror .eq. 1 .and. icontorno .lt. n-3 )
c	  icontorno=icontorno+3
c         call inicializa(icontorno,etal,s,ds,oo,sol,ierror)
c	end do
	  
	do i=icontorno,n
	   do jj=1,4
	      do ii=1,4
                ssum=0.
	         do k=1,4
	            ssum = ssum +  etal(ii,k,i)*etal(k,jj,i)
	         enddo
	         etal2(ii,jj,i) = ssum
	      enddo
	   enddo
	enddo

 	do i=icontorno+1,n
           delta=deltai(i)
           delt2=delt2i(i)
	   do jj=1,4
              e0=etal(jj,1,i-1)
	      e1=etal(jj,1,i)
              s00=s(i-1)
              s01=s(i)
	      den(jj) =-( e0*s00 + e1*s01 )*
     &	                 delta + sol(jj,i-1)
	      den(jj) =-(detal(jj,1,i-1)*s00 - detal(jj,1,i)*s01  +
     &                  e0*ds(i-1) -  e1*ds(i) +
     &                  etal2(jj,1,i-1)*s00 - etal2(jj,1,i)*s01 )*
     &                  delt2 + den(jj)

	      do ii=1,4
	         den2(ii,jj)=etal(ii,jj,i-1)*delta+
     & delt2*(detal(ii,jj,i-1)+etal2(ii,jj,i-1))
	         den(jj)=den(jj) + (etal(jj,ii,i-1)*delta+
     & delt2*(detal(jj,ii,i-1)+etal2(jj,ii,i-1)))*sol(ii,i-1)
	         num(ii,jj) = - etal(ii,jj,i)*delta + 
     & delt2*(detal(ii,jj,i)+etal2(ii,jj,i))
	      enddo
	      num(jj,jj) = 1.0 + num(jj,jj)
	      den2(jj,jj) = 1.0 + den2(jj,jj)
	   enddo

c	   do ll=1,4
c		write(*,*) 'hola hermite',(num(ll,lll),lll=1,4),delta,
c     & delt2,etal(ii,jj,icontorno+1),detal(ii,jj,icontorno+1),etal2(ii,jj,icontorno+1)
c	   enddo
             
        
	   call matinx(num)

	   do jj=1,4
              sum=0.
	      do ii=1,4
	         sum = sum + num(jj,ii)*den(ii)
	         sumoo= 0.e0
	         do k=1,4
	            sumoo = sumoo + num(ii,k)*den2(k,jj)
	         enddo
	         OO(ii,jj,i-1) = sumoo
	      enddo
	      sol(jj,i) = sum
	   enddo
	enddo

	do ii=1,4
	   svec(ii) = sol(ii,n)
	enddo

	do jj=1,4
	   do ii=1,4
	      OO(ii,jj,n) = 0.e0
	   enddo
	   OO(jj,jj,n) = 1.e0
	enddo

	do i=n-2,1,-1
	   do jj=1,4
	      do kk=1,4
                 sum=0.e0
	         do k=1,4
	            sum=sum+OO(kk,k,i+1)*OO(k,jj,i)
	         enddo
	         den2(kk,jj) = sum
              enddo
	   enddo
	   do jj=1,4
	      do kk=1,4
	         OO(kk,jj,i)=den2(kk,jj)
	      enddo
	   enddo
	enddo

	do i=1,n
	   s0(i) = 0.e0
	   sol(1,i) = sol(1,i)-s(i)
	enddo

	if(mnodos(1).ne.0)call DELTA1(n,tk,bt,sol,OO,etall,rt)
c	if(mnodos(2).ne.0)call DELTA1(n,pk,s0,sol,OO,etall,rp)
	if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &     mnodos(6).ne.0)call DELTA1(n,pk,s0,sol,OO,etall,rp)
	if(mnodos(3).ne.0)call DELTA1(n,mk,s0,sol,OO,etall,rm)
	if(mnodos(4).ne.0)call DELTA1(n,hk,s0,sol,OO,etall,rh)
	if(mnodos(5).ne.0)call DELTA1(n,vk,s0,sol,OO,etall,rv)
	if(mnodos(6).ne.0)call DELTA1(n,gk,s0,sol,OO,etall,rg)
	if(mnodos(7).ne.0)call DELTA1(n,fk,s0,sol,OO,etall,rf)

	return
	end

c ___________________________________________________________________
       subroutine inicializa(icontorno,etal,s,ds,oo,sol,ierror)

       implicit real*4 (a-h,o-z)

       include 'PARAMETER' !por kt
       parameter (kt16=16*kt)
       real*4 etal(4,4,kt),inveta(4,4),s(*),ds(*),oo(4,4,kt),sol(4,kt)

       do k=1,icontorno
	  do j=1,4
	     do i=1,4
	        inveta(i,j)=etal(i,j,k)
c		 print*,'inveta(i,j)=',i,j,inveta(i,j)
	     end do
	  end do

	  call matinx2(inveta,ierror)	!calcula la matriz inversa 4*4

	  do i=1,4
c	     print*,'en hermite.f, inveta(i,k)=',i,k,inveta(i,k),ds(k)
	     sol(i,k)=inveta(i,1)*ds(k)
	  end do
	  sol(1,k)=sol(1,k)+s(k)
       end do

	do l=2,icontorno+1   !aprox. de difusion para el op. evolucion
	   do j=1,4
              do i=1,4
                OO(i,j,l)=0.
	      enddo
c	      if(i.eq.j)OO(i,j,l)=0. !1.
            enddo       
c           smed=.5*(s(l)-s(l-1))
c	    smed=s(l)
	    OO(1,1,l-1)=0.   !(sol(1,l)-smed)/(sol(1,i-1)-smed)
	enddo

       return
       end
