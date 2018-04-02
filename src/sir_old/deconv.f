c DECONV convoluciona vin con una gaussiana de "velocidad" s
c y con el filtro leido del fichero filtro (col. 1 ldo, col2 trans.)
c da el resultado en vin
c isigno (=1 convolucion),(=-1 deconvolucion),(0= return)
c n es el numero de puntos
c el vector a deconvolucionar entra por medio del common vobservado

	subroutine deconv(vin,isigno,nlins,npass,dlamda0s,dlamdas,s)

	parameter(m1=1024)   !m1 must be an integer power of 2!!!! 

	real*4 vin(*),dlamda0s(*),dlamdas(*)
	real*4 frec,ex,cota,paso,sigma,pi,c,s
	real v1(m1),v2(m1),expo(m1)
	integer npass(*),ifiltro
	common/ifiltro/ifiltro
	data ivez/0/

        ntot=m1  !antes ponia ntot=kld

	c=3.e5		!c en km/seg luego "s vendra en km/seg"
	pi=3.141590    !dlamda0 vendra dada en a
	cota=26.
        
	if(ivez.eq.0.and.ifiltro.eq.1)then
          ivez=1
          call psf(ntot,nlins,npass,dlamda0s,dlamdas) !da la transf. en ftrans
	end if

	k=0
	if(isigno.eq.0.and.ifiltro.ne.1)return
	if(s.eq.0.and.ifiltro.ne.1)return

c	print*,'deconv: nlins=',nlins,' vmac=',s,' ntot=',ntot

c descomponemos vin en datos para cada linea: v1
	k1=0
	kk=0
   
c rellenamos hasta 'ntot' con ceros
	do j=1,nlins
	  n2=npass(j)
	  i1=(ntot-n2)/2
	  i2=i1+n2

          if(n2.gt.ntot)then
            print*,'en DECONV: la linea ',j,' tiene mas de ',ntot,' frecuencias'
            stop
          end if

	  if(n2.gt.1)then
	    sigma=1.e3*dlamda0s(j)*s/c
	    paso=abs((dlamdas(k1+2)-dlamdas(k1+1)))
c	    print*,'deconv: npass(',j,')=',npass(j),' sigma=',sigma,' paso=',paso
c	   print*,'_____: dlamda0s(',j,')=',dlamda0s(j),' sigma=',sigma,' paso=',paso
	    if(paso.lt.1.e-20)paso=1.
	    do i=1,npass(j)
		k1=k1+1
		v1(i)=vin(k1)		!*1.e-13
	    end do

	    do i=1,i1
	       v2(i)=v1(1)
	    end do

	    do i=i2+1,ntot
	       v2(i)=v1(n2)
	    end do

	    do i3=i1+1,i2
		v2(i3)=v1(i3-i1)
	    end do

	    do i=1,ntot/2
		frec=(i-1)/(paso*ntot)
	        ex=2.*pi*pi*sigma*sigma*frec*frec
		if(ex.gt.cota)ex=cota
	        expo(i)=exp(-1.*isigno*ex)
	    end do
	    do i=ntot/2+1,ntot
		frec=(ntot-(i-1))/(paso*ntot)
	        ex=2.*pi*pi*sigma*sigma*frec*frec
		if(ex.gt.cota)ex=cota
	        expo(i)=exp(-1.*isigno*ex)
	    end do

	    linea=j 
	    if(n2.gt.1)call deconvsub2(v2,ntot,expo,linea)

	    do i=i1+1,i2
		kk=kk+1
		vin(kk)=v2(i)		
	    end do
	  else
	    do i=1,n2
	       kk=kk+1
	       k1=k1+1
	    end do
	  end if
	end do

	return
	end
c_______________________________________________________________________
c psf lee la PSF del fichero filtro se supone que contiene dos columnas
c la primera con la ldo en mA centrada en el maximo
c la segunda la psf (no necesita estar normalizada en area)
c escribe su transformada en ftrans 
c da los valores solo para el muestreo definido por dlamda para  
c cada linea 
c num=numero de puntos del fichero filtro
c ntot=numero de puntos que quiero

	subroutine psf(ntot,nlin,npas,dlamda0,dlamda)

	include 'PARAMETER'  !solo por kl

	parameter (m1=1024,m2=2*m1,nmx=401)
	real*4 dlamda0(*),dlamda(*)
        real*4 x(nmx), y(nmx)
	
	real*4 entrada(m2)
	real*4 ftrans(4*kl,m2)
	real*4 xa(11),ya(11)
	integer npas(*)

        character*100 filtro,control
    
        common/filtro/filtro,x,y,num
	common/canal/icanal
	common/nombrecontrol/control
	common/ftransformada/ftrans

        ngrado=2        !interpolo con parabolas
	n2=int(ngrado/2)


c leemos el fichero filtro
c 	open(55,file=filtro,status='old',err=100)
c 	do ii=1,nmx
c            read(55,*,end=127,err=101)x(ii),y(ii)
c 	end do
c 127	num=ii-1

	if(num.eq.nmx-1)print*,'The psf file is being truncated (it has more than 401 wavelengths!)'


c interpolamos
	k1=0 
        kcount=0  
        ntot2=ntot/2
	do j=1,nlin
	   n=npas(j)
	   paso=abs((dlamda(k1+2)-dlamda(k1+1)))
	   k1=k1+n
	  	   
	   do i=1,ntot
c	      x0=(i-1-ntot/2)*paso
              if(i.le.ntot2)then
                 x0=(i-1)*paso
              else
                 x0=(i-1-ntot)*paso
              end if
        
              entrada(2*i)=0.

	      if(x0.lt.x(1).or.x0.gt.x(num))then
		 entrada(2*i-1)=0.
              else
	         call locate(x,num,x0,jj)  !jj es el indice anterior
	         n3=jj-n2-1
		 if(n3.lt.0)n3=0
                 if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	         do k=1,ngrado+1
		    xa(k)=x(n3+k)
	         end do
	         do k=1,ngrado+1
		    ya(k)=y(n3+k)
	         end do
	         call polint(xa,ya,ngrado+1,x0,salida,error)
                 entrada(2*i-1)=salida
	      end if
c              print*,i,entrada(2*i-1)
           end do
	   
 	   call four1(entrada,ntot,1)

           area=entrada(1)
	   do i=1,2*ntot
	      ftrans(j,i)=entrada(i)/area
	   end do 
 
	end do !do en lineas
        close(55)      
	

	return

100	print*,' '
	print*,'STOP: The file containing the PSF does NOT exist:'
        print*,filtro
        print*,'_______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*),'STOP: The file containing the PSF does NOT exist:'
c        write(icanal,*) filtro
c        close(icanal)
	stop

101     print*,' '
	print*,'STOP: The file containing the PSF is not written correctly:'
        print*,filtro
        print*,'_______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*),'STOP: The file containing the PSF is not written correctly:'
c        write(icanal,*) filtro       
c        close(icanal)
	stop
 

	end
c_______________________________________________________________________




c_______________________________________________________________________
c DECONVSUB2 (incluida en DECONV) multiplica las transformadas de
c Fourier de la gaussiana, de la PSF y del perfil y anti-transforma

c rutina deconvsub2(data,n,datb,isigno)
c isigno=1 convolucion, -1 deconvolucion
c data : vector a convolucionar
c datb : transformada del vector con/deconvolucionador

	subroutine deconvsub2(data,n,datb,linea)

	include 'PARAMETER'
	parameter(m1=1024,m2=2*m1)

	real*4 data(*),datb(*),entrada(m2)
	real*4 ftrans(4*kl,m2)
	integer ifiltro
	common/ftransformada/ftrans
	common/ifiltro/ifiltro

	if(n.le.1)return
	if(n.ne.m1) stop 'n not equal to m1 in deconvsub2'
        
        do i=1,n
           entrada(2*i-1)=data(i)
           entrada(2*i)=0.
        end do
 	call four1(entrada,n,1)
	
	k=0
	if (ifiltro.eq.1)then
	   do i=1,n
              preal1=entrada(k+1)*datb(i)
              pimag1=entrada(k+2)*datb(i)
              preal2=ftrans(linea,k+1)
              pimag2=ftrans(linea,k+2)
              entrada(k+1)=preal1*preal2-pimag1*pimag2  !parte real
              entrada(k+2)=preal1*pimag2+preal2*pimag1  !parte imaginaria
              k=k+2
           end do 
	else
           do i=1,n
              entrada(k+1)=entrada(k+1)*datb(i) !parte real
              entrada(k+2)=entrada(k+2)*datb(i) !parte imaginaria
              k=k+2
           end do 
	end if
	   
	call four1(entrada,n,-1)

c normalizamos dividiendo por el numero de puntos
        xn=float(n) 

        do i=1,n
           data(i)=entrada(2*i-1)/xn !parte real
        end do 

	return
	end
	
c********************************************************
c	Subrutina FOUR1 del NUMERICAL RECIPES
	
	subroutine four1(data,nn,isigno)
	real*8 wr,wi,wpr,wpi,wtemp,theta
	dimension data(*)

c        print*,'entro en four1 00',nn,2*nn
	n=2*nn
	j=1
	do 11 i=1,n,2
c          print*,'estoy en four1 do 11',i

	  if (j.gt.i) then
	    tempr=data(j)
	    tempi=data(j+1)
	    data(j)=data(i)
	    data(j+1)=data(i+1)
	    data(i)=tempr
	    data(i+1)=tempi
	  endif
	  m=n/2
 1	  if ((m.ge.2).and.(j.gt.m)) then
	    j=j-m
	    m=m/2
	    goto 1
	  end if
	  j=j+m
 11	continue
c        print*,'estoy en four1 01' 

     

	mmax=2
 2	if (n.gt.mmax) then
c         print*,'estoy en four1 despues de goto 2'

	istep=2*mmax
	theta=6.28318530717959d0/(isigno*mmax)
	wpr=-2.d0*dsin(0.5d0*theta)**2
	wpi=dsin(theta)
	wr=1.d0
	wi=0.d0
	do 13 m=1,mmax,2
	  do 12 i=m,n,istep
	    j=i+mmax
	    tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
	    tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
	    data(j)=data(i)-tempr
	    data(j+1)=data(i+1)-tempi
	    data(i)=data(i)+tempr
	    data(i+1)=data(i+1)+tempi
 12	  continue
	  wtemp=wr
	  wr=wr*wpr-wi*wpi+wr
	  wi=wi*wpr+wtemp*wpi+wi
 13	continue
	mmax=istep
	goto 2
	endif
c        print*,'estoy en four1 99'

	return
	end
