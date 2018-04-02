c EQUISUBMU evalua la presion en equilibrio hidrostatico
c_______________________________________________________________
c equisubmu rutina que evalua la presion en equilibrio hidrostatico
c teniendo en cuenta el numero de electrones para el calculo de mu
c Luis Bellot y Basilio Ruiz 27/7/95 
c Cambio del metodo de integracion Basilio Ruiz 3/9/96
c _______________________________________________________________
	subroutine equisubmu(ntau,tau1,t,pe,pg,z,ro)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'   !solo por kt, que era igual a 1000
	parameter (nex=28,cgases=83145100.)
        real*4 tau1(kt)
	real*4 tau(kt),t(*),pe(*),pg(*),kac,d2,x(kt),kappa(kt),taue(kt)
        real*4 z1(kt),z(*),ro(*),y(kt)
	real*4 mu,pcontorno
	integer*4 nmaxitera
    
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth                   !esto esta deshabilitado (entra 1)
        common/precisoitera/precitera      
        common/anguloheliocent/mu
        common/nmaxitera/nmaxitera
c        common/pinicial/pcontorno
        
        precitera=1.e-5
        nmaxitera=250
        
	g=mu*2.7414e+4		!gravedad cm/s^2 en fotosfera solar   
       	avog=6.023e23
       	cth=1.
       	mu=1.

        do i=1,ntau
           tau(i)=tau1(i)+alog10(cth)  !this is the optical depth in 
        end do                        !the VERTICAL direction for 
	                              !hydrostatic equilib. computations
        tau0=1.
        imin=1
        do i=2,ntau
           if(abs(tau(i)).lt.tau0)then
              imin=i         !zero of z is at log tau (vertical) =0
              tau0=tau(i)
           end if
        end do

        do i=1,10
           d1(i)=0
        end do
c calculamos el peso molecular medio pmu
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do

	do i=1,ntau
           taue(i)=10.**tau(i)
        end do

	do i=1,ntau-1
            x(i+1)=g*(tau(i)-tau(i+1))*2.30259  
        end do

c inicializamos 
	tsi=t(ntau)
        psi=pe(ntau)
c        pg(ntau)=pcontorno
        call gasc(tsi,psi,pg(ntau),pp)
	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
        pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
        kappa(ntau)=kac*avog/pmusum
        ro(ntau)=pesomedio*pg(ntau)/tsi/cgases
        y(ntau)=taue(ntau)/kappa(ntau)/ro(ntau)

c integramos
        do i=ntau-1,1,-1
           pgpr=pg(i+1)+x(i+1)*taue(i+1)/kappa(i+1)
	   tsi=t(i)

           call pefrompg11(tsi,pgpr,psi)
           pe(i)=psi

           call gasc(tsi,psi,pgpr,pp)
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
           pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
           kappa(i)=kac*avog/pmusum

           pg(i)=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i))/2.

           nit=1
           dif=2.*abs((pg(i)-pgpr))/(pg(i)+pgpr)

           do while (nit.lt.20.and.dif.gt.1.e-5)
              nit=nit+1
              pgpr=pg(i)
              call pefrompg11(tsi,pgpr,psi)
              pe(i)=psi

              call gasc(tsi,psi,pgpr,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
              pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
              kappa(i)=kac*avog/pmusum
              pg(i)=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i))/2.
              dif=2.*abs((pg(i)-pgpr))/(pg(i)+pgpr)
            end do

            if(dif.gt.0.1)print*,'WARNING: 
     & Hydrostatic equilibrium results in inaccurate electron pressures '
            call pefrompg11(tsi,pg(i),psi)
            pe(i)=psi
            ro(i)=pesomedio*pg(i)/tsi/cgases
            y(i)=taue(i)/kappa(i)/ro(i)
	end do
        z1(1)=0.
        do i=2,ntau  !this is the z scale along the vertical direction
           z1(i)=z1(i-1)+x(i)/g*(y(i-1)+y(i))/2.
        end do 
        do i=1,ntau
           z(i)=(z1(i)-z1(imin))*1.e-5
        end do

c        do i=1,ntau   !we output geometrical heights ALONG the LOS
c          z(i)=z(i)/mu
c        enddo


       return
       end 

c____________________________________________________________________
c PEFROMPG10 (incuida en EQUISUBMU) evalua la presion electonica  desde la 
c temperatura y la presion gaseosa 
c pefrompg10 evalua la presion electonica p correspondiente a t1 y pg
c
c	subroutine pefrompg10(t,pg,p)
c
c	implicit real*4 (a-h,o-z)
c        common/preciso/prec      
c 
c        dif=1.
c        n2=0
c        p1=p
c        do while (dif.gt.prec.and.n2.lt.50)
c	     n2=n2+1
c	     p=(p+p1)/2.
c             p1=p
c             call pe_pg10(t,p,pg)
c             dif=abs(abs( (p-p1)/(p+p1) )-1.0)
c         end do
c	return
c	end
c____________________________________________________________________
c PEFROMPG11 (incuida en EQUISUBMU) evalua la presion electonica  desde la 
c temperatura y la presion gaseosa 
c pefrompg11 evalua la presion electonica p correspondiente a t1 y pg

	subroutine pefrompg11(t,pg,p)

	implicit real*4 (a-h,o-z)
        common/precisoitera/prec   
        common/nmaxitera/nmaxitera
 
        dif=1.
        n2=0
        
         p0=p
         call pe_pg10(t,p,pg)
         dif=abs( (p-p0)/p )
         if(dif .gt. 1)call inicia_pefrompgt(t,pg,p)

        p1=p
        do while (dif.gt.prec.and.n2.lt.nmaxitera)
	     n2=n2+1
	     p=(p+p1)/2.
             p1=p
             call pe_pg10(t,p,pg)
             dif=abs( (p-p1)/p )
         end do
c         print*,'pefrompg11 n2=',n2,dif,p,p1
	return
	end
c____________________________________________________________________
c PE_PG10 (incuida en EQUISUBMU) evalua la presion electonica  desde la 
c la presion gaseosa y una estimacion de la presion electronica 
c calcula la presion electronica a partir de la pg y de una estimacion de la pe

      subroutine pe_pg10(t,pe,pg)

      parameter (ncontr=28)
      real*4 cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)

      real*4 du0,du1,du2,dcmol(91)

    
      if(t.lt.500)then
         print*,'pe_pg10: temperature < 500 K '
	 print*,'temperature = 500 K'
	 t=500.
      end if	 
      theta=5040./t
      g4=0.
      g5=0.
      if(pe.le.0)then
         pe=1.e-15
         g4=0.
         g5=0.
      else
         call molecb(theta,cmol,dcmol)
         do i=1,2
            call acota(cmol(i),-30.,30.)
         end do
         g4=pe*10.**(cmol(1))
         g5=pe*10.**(cmol(2))
      end if 
c ahora calculo los niveles u0,u1,u2 y sus derivadas
      do 5 i=1,ncontr
         iii=i
5     	 call neldatb(iii,0.,weight,alfai(i),chi1(i),chi2(i))
6     do 4 i=1,ncontr
      	 iii=i
4     	 call nelfctb(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2)

 
      g2=saha(theta,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
      
      g3=saha(theta,0.754,1.,u0(1),pe)        ! p(h)/p(h-) 

      call acota(g3,1.e-30,1.e30)
      g3=1.d0/g3                              ! p(h-)/p(h) 
     
      g1=0.
      do 1 i=2,ncontr
        a=saha(theta,chi1(i),u0(i),u1(i),pe)
        b=saha(theta,chi2(i),u1(i),u2(i),pe)
	c=1.+a*(1.+b)
1	g1=g1+alfai(i)/c*a*(1.+2.*b)
      
        a=1.+g2+g3
        b=2.*(1.+g2/g5*g4)
        c=g5
        d=g2-g3
        e=g2/g5*g4

	call acotasig(a,1.e-15,1.e15)
	call acotasig(d,1.e-15,1.e15)

        c1=c*b**2+a*d*b-e*a**2

        c2=2.*a*e-d*b+a*b*g1
	
        c3=-(e+b*g1)

        f1=0.5*c2/c1

        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1)
        f5=(1.-a*f1)/b
        f4=e*f5
        f3=g3*f1
        f2=g2*f1
        fe=f2-f3+f4+g1

        call acota(fe,1.e-30,1.e30)
	phtot=pe/fe

        if(f5.gt.1.e-4) goto 2
          const6=g5/pe*f1**2
          const7=f2-f3+g1
	  do 3 i=1,5
      	       f5=phtot*const6
	       f4=e*f5
	       fe=const7+f4
	       phtot=pe/fe
3	  continue
      
2       pe=pg/(1.+(f1+f2+f3+f4+f5+0.1014)/fe)
        if(pe.le.0)pe=1.e-15

      return
      end
c ____________________________________________________________________________
      
c esta rutina interpola en ttau el modelo tau,t,pe y da la salida en
c ttau,tt,ppe por medio de un polinomio de grado ngrado

	subroutine interpolatp(ngrado,n,tau,t,pe,nnew,ttau,tt,ppe)

	real tau(*),t(*),pe(*),ttau(*),tt(*),ppe(*)
	real xa(11),ya(11)

c interpolaremos las presiones en logaritmos neperianos
	num=n
	do i=1,num
	   pe(i)=alog(pe(i))
	end do

c interpolamos
	n2=int(ngrado/2)
	
	do i=1,nnew
	   CALL LOCATE(TAU,NUM,TTAU(I),J)
	   n3=j-n2-1
           if(n3.lt.0)n3=0
           if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	   do k=1,ngrado+1
	      xa(k)=tau(n3+k)
	   end do
	   do k=1,ngrado+1
	      ya(k)=t(n3+k)
	   end do
	   CALL POLINT(XA,YA,NGRADO+1,TTAU(I),TT(I),ERROR)

	   do k=1,ngrado+1
	      ya(k)=pe(n3+k)
	   end do
	   CALL POLINT(XA,YA,NGRADO+1,TTAU(I),ppe(I),ERROR)
        end do

        
	do i=1,nnew
           ppe(i)=exp(ppe(i))
        end do
        do i=1,n
           pe(i)=exp(pe(i))
	end do

	return
	end
c ____________________________________________________________________________

	subroutine acota(x,x0,x1)

	if(x.lt.x0)x=x0	
        if(x.gt.x1)x=x1
   
        return
        end 
c ____________________________________________________________________________
       
	subroutine acotasig(x,x0,x1)

        if(x.lt.0)then
           x=-x
           call acota(x,x0,x1)
           x=-x
         else
           call acota(x,x0,x1)
         end if
        return
        end 
c ____________________________________________________________________________
     