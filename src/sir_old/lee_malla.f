c rutina lee_malla
c lee el fichero de control de lineas y longitudes de onda:'mallaobs'

c ntl   : numero total de lineas
c nli   : numero total de frecuencias ,incluyendo todas las lineas
c nlin  : indice de cada linea y de cada blend
c npas  : numero de puntos en cada linea
c dlamda: cada delta de l.d.o. en mA
c nble  : numero de blends de cada linea

	subroutine lee_malla(mallaobs)
	
	include 'PARAMETER'  !por kl
	parameter (kl4=4*kl,kld4=4*kld)
	
	integer ntl,nlin(kl),npas(kl),nble(kl)
	real*4 dlamda(kld)
	integer ntls,nlins(kl4),npass(kl4)
	real*4 dlamdas(kld4)       
	integer nd(kl)
        integer blanco,ifail,ierror,ntonto,jj,k
	real*4 dini,dipa,difi,errpasos,numeropasos
	character*100 mallaobs
        character*31 mensajito
	character*80 men1,men2,men3
     
        common/Malla/ntl,nlin,npas,nble,dlamda  !common para StokesFRsub
        common/Malla4/ntls,nlins,npass,dlamdas  !common para StokesFRsub

	men1=' '
	men2=' '
	men3=' '
       
	ican=57
	mensajito=' containing the wavelength grid'
	call cabecera(ican,mallaobs,mensajito,ifail)
	if(ifail.eq.1)goto 999

        jj=0
        k=0 
        nli=0
	numlin=0
	do i=1,1000   !Vamos a tener mas de 1000 lineas???? Pues entonces....  
	  
	  call mreadmalla(ican,ntonto,nd,dini,dipa,difi,ierror,blanco)
	  if(ierror.eq.1)goto 8
	  if(blanco.ne.1)then

	    numlin=numlin+1   
	    if(numlin.gt.kl.or.jj+ntonto.gt.kl)then !comprobamos numero de lineas
	     men1='STOP: The number of lines in the grid is larger than the '
	     men2='limit. Decrease this number or change the PARAMETER file.'
	     call mensaje(2,men1,men2,men3)
            end if
	    nble(numlin)=ntonto

	    do j=1,nble(numlin)
	      jj=jj+1
              nlin(jj)=nd(j)
	    end do

	    numeropasos=(difi-dini)/dipa+1
	    npas(numlin)=int(numeropasos)
	    errpasos=10*(numeropasos-int(numeropasos))	    
            if(errpasos.gt.0.5) npas(numlin)=npas(numlin)+1
            if(nli+npas(numlin).gt.kld)then  !comprobamos numero longitudes onda
	       men1='STOP: The grid has more points than the limit kld.'
	       men2='Decrease the wave number or change the PARAMETER file.'
	       call mensaje(2,men1,men2,men3)
            end if
	    do j=1,npas(numlin)
	       nli=nli+1
	       dlamda(nli)=dini+(j-1)*dipa
	    end do
         endif
	enddo
8	ntl=numlin

c	print*,'Number of wavelengths in the wavelength grid: ',nli
	close(ican)

c	print*,'lee_malla',ntl,nlin(ntl),npas(ntl),nble(ntl),dlamda(1)
	
	
        ntls=0
        nb=0
        k4c=0
	do i=1,4
	   jj=0
	   k3c=0
	   do k=1,ntl
              do jble=1,nble(k)   
                 nb=nb+1
                 jj=jj+1
                 nlins(nb)=nlin(jj)
	       end do
               ntls=ntls+1
	       npass(ntls)=npas(k)
	       do l=1,npas(k)
	           k3c=k3c+1
	           k4c=k4c+1
	           dlamdas(k4c)=dlamda(k3c)
	       end do 
	    end do   
	end do
	

	return

c       -------------------------------------------------------------------------
c	Mensajes de error:

999	men1=' '          !'STOP: The grid file does NOT exist:'
	men2=mallaobs
	call mensaje(2,men1,men2,men3)

992	men1=' '          !'STOP: Incorrect format in grid file:'
        men2=mallaobs
	call mensaje(2,men1,men2,men3)


	end
