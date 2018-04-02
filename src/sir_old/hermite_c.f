	subroutine hermite_c(s,ds,etall,n,svec,nn,bt,tk,pk,
     &           vk,mk,rt,rp,rv,rm)

        implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kt
        parameter (kt16=16*kt)
        parameter (dl10=2.3025851)

	integer n,nn
	real*4 svec

c internas
        real*4 sol(kt),tau(kt),taue(kt)
        real*4 deltae(kt),deltai(kt),delt2i(kt)
	real*4 ds(*), detal(kt), etal(kt), etall(*)
	real*4 etal2(kt)
	real*4 s(*), den, num, den2
	real*4 OO(kt)
        real*4 tk(kt),pk(kt),vk(kt),mk(kt)
	real*4 bt(kt)
	real*4 rt(*),rp(*),rv(*),rm(*)
c       deltae(i)=(taue(i)-taue(i+1))
c       deltai(i) = (tau(i)-tau(i-1))/2.
c       delt2i(i) = delta(i)*delta(i)/3. 
        common/segunda/tau,taue,deltae,deltai,delt2i

c se calcula donde poner el contorno
        icontorno=1
	do i=1,n-1
           taunu=abs(etall(i))*deltae(i)
           if(taunu.gt.5.)icontorno=i
        end do
c	print*,'en hermite_c icontorno=',icontorno
c	icontorno=1
c se calcula la derivada de la matriz de absorcion
	do i=1,n
	   etal(i) = etall(i)*taue(i)*dl10
  	end do 	
	call derivacuad(etal,detal,n,tau)

c se integra
        call inicializa_c(icontorno,etal,s,ds,oo,sol,ierror)
	if(ierror.eq.1)return
	do i=icontorno,n
  	   etal2(i) = etal(i)*etal(i)
	enddo

c	sol(1) = s(1) + ds(1)/etal(1)

 	do i=icontorno+1,n
           delta=deltai(i)
           delt2=delt2i(i)
           e0=etal(i-1)
	   e1=etal(i)
           s00=s(i-1)
           s01=s(i)
	   den =-( e0*s00 + e1*s01 )*
     &	                 delta + sol(i-1)
	   den =-(detal(i-1)*s00 - detal(i)*s01  +
     &                  e0*ds(i-1) -  e1*ds(i) +
     &                  etal2(i-1)*s00 - etal2(i)*s01 )*
     &                  delt2 + den

	   den2=e0*delta+
     & delt2*(detal(i-1)+etal2(i-1))
	   den=den + (e0*delta+
     & delt2*(detal(i-1)+etal2(i-1)))*sol(i-1)
	   num= - e1*delta +
     & delt2*(detal(i)+etal2(i))
	   num = 1.0 + num
	   den2 = 1.0 + den2
        
	   OO(i-1) = den2/num
	   sol(i) = den/num
	enddo

	svec = sol(n)
	OO(n) = 1.0

	do i=n-2,1,-1
	   OO(i)=OO(i+1)*OO(i)
	enddo

c calculamos las f.resp
	do i=1,n-1
           ss=sol(i) - s(i)
           ooo=OO(i)

           rt(i)=ooo*(tk(i)*ss-etall(i)*bt(i))
           rp(i)=ooo*pk(i)*ss
           rv(i)=ooo*vk(i)*ss
           rm(i)=ooo*mk(i)*ss
c           print*,'hermite_c sol,ooo,tk,ss,etall,bt=',sol(i),ooo,tk(i),ss,etall(i),bt(i)
	enddo
c	print*,'hermite_c sol(n)=',svec
	rt(n)=0.
	rp(n)=0.
	rv(n)=0.
	rm(n)=0.

	return
	end
c ________________________________________________________________________
       subroutine inicializa_c(icontorno,etal,s,ds,oo,sol,ierror)

       implicit real*4 (a-h,o-z)

       include 'PARAMETER' !por kt
       parameter (kt16=16*kt)
       real*4 etal(*),inveta,s(*),ds(*),oo(*),sol(*)

       do k=1,icontorno
	 inveta=0.
	 if(abs(etal(k)) .gt. 1.e-30)then inveta=1./etal(k)
	 sol(k)=inveta*ds(k)+s(k)
	 OO(k)=0.
       end do

       return
       end

