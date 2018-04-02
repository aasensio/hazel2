c taulinea
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau : numero de puntos en tau
c tau  : log10(tau)
	subroutine taulinea(il,cth,ihe,vx,atmos,ntau)

	include 'PARAMETER'  !por kt
	implicit real*4 (a-h,o-z)
	real*4 atmos(*),sqrtb
	real tau(kt),t(kt),p(kt),vmic(kt),h(kt),vz(kt),gamma(kt),phi(kt)
	real ttau(kt)
	integer il,ihe

	pi=3.14159265
	pi180=pi/180.
	epsilon=1.e-3

              if(il.eq.0)then
		do i=1,ntau
	          atmos(i+6*ntau)=atmos(i+6*ntau)*pi180
	          atmos(i+7*ntau)=atmos(i+7*ntau)*pi180
                enddo
c		print*,'il taulinea gamma y fi',atmos(6*ntau+1),atmos(7*ntau+1),il
             else
		do i=1,ntau
	          atmos(i+6*ntau)=atmos(i+6*ntau)/pi180
	          atmos(i+7*ntau)=atmos(i+7*ntau)/pi180
c		  print*,'il taulinea gamma y fi',atmos(6*ntau+1),atmos(7*ntau+1),il
                enddo
             endif



        if(cth.eq.1.) return



	if(abs(cth).ge.1)then
	   sth=0.
	else
	   sth=ihe*sqrtb(1.e0-cth*cth)	!seno de teta
	end if

	do i=1,ntau
	   ttau(i)=atmos(i)
	   t(i)=atmos(i+ntau)
	   p(i)=atmos(i+2*ntau)
	   vmic(i)=atmos(i+3*ntau)
	   h(i)=atmos(i+4*ntau)
	   vz(i)=atmos(i+5*ntau)
	   gamma(i)=atmos(i+6*ntau)
	   phi(i)=atmos(i+7*ntau)
	end do

	if(il.eq.0)then
	   do i=1,ntau
	      tau(i)=ttau(i)-alog10(cth)
              vz(i)=vx*sth+vz(i)*cth
	      cgaz=cos(gamma(i)*pi/180.)
	      sgaz=sin(gamma(i)*pi/180.)
	      cfiz=cos(phi(i)*pi/180.)
	      sfiz=sin(phi(i)*pi/180.)
	      cga=-sth*cfiz*sgaz+cth*cgaz
	      if(abs(cga).ge.1.)then
                  print*,'coseno de gamma mayor que uno'
		  cga=1.-epsilon/100.
	      endif
	      sga=sqrtb(1.e0-cga*cga)

	      if(sgaz.lt.0)sga=-sga
	      if(abs(sga).le.epsilon)then
	         gamma(i)=0.    !gamma
	         phi(i)=0.             !atencion, cambio para HORST
	      else
	         gamma(i)=atan2(sga,cga)
c		 gamma(i)=acos(cga)

	         sfi=sfiz*sgaz/sga
	         cfi=(cth*cfiz*sgaz+sth*cgaz)/sga
	         phi(i)=atan2(sfi,cfi)
              end if

              ggg=gamma(i)
	      if(ggg.lt.-2.*pi)ggg=ggg-2.*pi*int(ggg/(2.*pi))
	      if(ggg.gt.-2.*pi.and.ggg.lt.-pi)ggg=2.*pi+ggg
	      if(ggg.lt.0.and.ggg.gt.-pi)ggg=-ggg
	      if(ggg.ge.2.*pi)ggg=ggg-2.*pi*int(ggg/(2.*pi))
	      if(ggg.gt.pi)ggg=2.*pi-ggg
              gamma(i)=ggg

	    end do

	    call quitasaltos2(ntau,pi,phi)

	  else
            do i=1,ntau
	       tau(i)=ttau(i)+alog10(cth)
	       vz(i)=(vz(i)-vx*sth)/cth
	       cga=cos(gamma(i))

	       sga=sin(gamma(i))
	       cfi=cos(phi(i))
	       sfi=sin(phi(i))
	       cgaz=sth*cfi*sga+cth*cga

	       if(abs(cgaz).ge.1.)cgaz=1.-epsilon/100.
	       sgaz=sqrt(1.e0-cgaz*cgaz)

	       if(sga.lt.0)sgaz=-sgaz
	       if(abs(sgaz).le.epsilon)then
	          gamma(i)=0.  !gamma
	          phi(i)=0.  !fi               !HORST
	         phi(i)=atmos(i+7*ntau)   !fi !atencion, cambio para HORST
	       else
	          gamma(i)=atan2(sgaz,cgaz)*180./pi
	          sfiz=sfi*sga/sgaz
	          cfiz=(cth*cfi*sga-sth*cga)/sgaz
	          phi(i)=atan2(sfiz,cfiz)*180./pi
               end if

                  ggg=gamma(i)
	          if(ggg.lt.-360.)ggg=ggg-360.*int(ggg/360.)
	          if(ggg.gt.-360..and.ggg.lt.-180.)ggg=360.+ggg
	          if(ggg.lt.0.and.ggg.gt.-180.)ggg=-ggg

 	          if(ggg.ge.360.)ggg=ggg-360.*int(ggg/360.)
	          if(ggg.gt.180.)ggg=360.-ggg
                  gamma(i)=ggg

	    end do

	    call quitasaltos2(ntau,180.,phi)
	
	  end if


c         interpolamos
	  call intmodel(3,ntau,ttau,tau,t,p,vmic,h,vz,gamma,phi)

	    do i=1,ntau
	       atmos(i)=tau(i)
	       atmos(i+ntau)=t(i)
	       atmos(i+2*ntau)=p(i)
	       atmos(i+3*ntau)=vmic(i)
	       atmos(i+4*ntau)=h(i)
	       atmos(i+5*ntau)=vz(i)
 	       atmos(i+6*ntau)=gamma(i)
	       atmos(i+7*ntau)=phi(i)
	    end do
	return
	end
c_____________________________________________________________________
c sqrtb funcion que calcula la raiz cuadrada
	real*4 function sqrtb(x)
	implicit real*4 (a-h,o-z)

	if(x.gt.0.e0)then
	   sqrtb=sqrt(x)
	else if(x.lt.-1.e-2)then
	   stop "raiz de algo menor que 0"
	else
	   sqrtb=0.e0
        end if
	return
	end
c____________________________________________________________________
c Rutina que quita los saltos en phi
	subroutine quitasaltos(ntau,cantidad,phi)

	real cantidad, phi(*)
	integer ntau

c	Se supone que phi esta ya en el rango [0,2pi]
c	Cantidad va a ser 180 grados o pi radiantes, dependiendo
c	de il

	    do i=1,ntau-1
		if (phi(i+1).lt.0..or.phi(i+1).gt.2.*cantidad)then 
 		   print*,'hay uno que se sale',i+1,phi(i+1)
		endif
		diff=phi(i+1)-phi(i)
c                if (abs(diff).gt..98*cantidad.and.abs(diff).lt.2.*.98*cantidad)then 
                if (abs(diff).gt..98*cantidad.and.abs(diff).lt.1.02*cantidad)then 
	              phi(i+1)=phi(i+1)-cantidad*diff/abs(diff)
		else 
                   if(abs(diff).ge.2.*.98*cantidad)then
	              phi(i+1)=phi(i+1)-2.*cantidad*diff/abs(diff)
		   endif
		endif		      			
	    enddo
	    icheck=0
	    do i=1,ntau
		if(phi(i).gt.2.*cantidad)icheck=1   
		if(phi(i).lt.0.)icheck=2 
	    enddo

	    isig=0
	    if(icheck.eq.1)then
		isig=-1
	    else if(icheck.eq.2)then
		isig=1
	    endif

	    do i=1,ntau
		phi(i)=phi(i)+isig*cantidad
	    enddo

	return
	end


c____________________________________________________________________
	subroutine quitasaltos2(ntau,cantidad,phi)


	real cantidad, phi(*)
	integer ntau

	do i=2,ntau
 	   diff=phi(i)-phi(i-1)  !esto tiene que ser un multiplo de cantidad
	   if(diff.eq.0.)then
		xquito=0.
	   else
		xquito=int(diff/(.95*cantidad))*diff/abs(diff)*cantidad
	   endif

           phi(i)=phi(i)+xquito
	   if(xquito.ne.0.)print*,'i,quito=',i,xquito
	enddo    

c       Se supone que ahora tenemos todos los fis sin saltos, referidos al primer
c	punto. Hay que desplazar la estratificacion para que el primer punto este
c	entre 0 y 2pi


	xmax=-10000.
	xmin=10000.

	do i=1,ntau
		if (phi(i).gt.xmax)xmax=phi(i)
		if (phi(i).lt.xmin)xmin=phi(i)
	enddo
c	print*,'xmax,xmin=',xmax,xmin

	return


	end


           

