c lee_model

	subroutine lee_model(modelname)

	include 'PARAMETER'  !por kt	
        parameter (kt8=8*kt+2)
        
        real*4 voffset,xmu
	integer ican,ntau
        integer*4 mnodos(18)	
	real*4 atmosmodel(kt8),a(8),pesostray
        real*4 tau(kt),t(kt),pe(kt),pg(kt),z(kt),ro(kt)
        character*100 modelname
	
	common/Atmosmodel/atmosmodel,ntau !common para StokesFRsub
	common/numeronodos/mnodos         !para StokesFRsub
        common/offset/voffset             !para StokesFRsub
        common/anguloheliocent/xmu        !para StokesFRsub
        
	ican=52

c offset de velocidad para perturbaciones relativas necesitamos que la velocidad sea siempre positiva        
        voffset=-15.e5    !cm/s
        xmu=1.            !coseno del angulo heliocentrico	
	
c contamos las lineas del modelo 1
	open(ican,file=modelname,status='old',err=991)
	ntau=0
	read(ican,*,err=999)a(1),fill1,peso1
	do while(ntau.lt.kt)
	    ntau=ntau+1
	    read(ican,*,end=221,err=999)a
	end do
221	ntau=ntau-1
        close(ican)
c       ahora leemos el modelo 
	open(ican,file=modelname)
	read(ican,*,err=999)atmosmodel(8*ntau+1),atmosmodel(8*ntau+2),pesostray
	do i=1,ntau
	   read(ican,*,err=999)(atmosmodel(i+j*ntau),j=0,7)
	end do
	close(ican)

c pasamos los angulos a radianes
        call taulinea(0,1.,1,0.,atmosmodel,ntau)
	
c definimos los nodos en todos los puntos (excepto para ls presion elctronica)
	do i=1,8                 
          mnodos(i)=ntau
        end do  
        mnodos(2)=0  

        call equisubmu(ntau,tau,t,pe,pg,z,ro)
 
        do i=1,ntau
            atmosmodel(i+2*ntau)=pe(i)
        end do
          
        return
        
991     print*,'STOP: The file containing model does NOT exist.'
	stop
        return
        
999     print*,'STOP: Incorrect format in the model file '
	stop
	return
	
	end
