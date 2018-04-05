c ________________________________________________________________
c rutina leelineasii
c datos de la linea

c Se supone que el fichero nomlineas tiene una cabecera con DOS
c lineas de menos de 80 caracteres

c       iln es el numero de la linea dentro del fichero lineas
c       atom es el simbolo del atomo,istage el estado de ionizacion,
c       wlengt la l.d.o. en angstroms
c       zeff es la correccion multiplicativa al damping (gamma6) 
c       energy potencial de excitacion del nivel inferior en ev.
c       loggf=log10(peso estadistico (g)* fuerza oscilador (f))
c       mult(i) da el spin total (spin=.5*float(mult(i)-1)),
c       design(i) da el momento angular orbital (ej:S=0,P=1,D=2,F=3...),
c       tam(i) da el momento angular total.
c       code1 representa mediante letras minusculas los momentos angulares
c       orbitales semienteros, asi p=1/2, f=3/2, h=5/2,k=7/2,m=9/2,o=11/2
c                              r=13/2,t=15/2, u=17/2,v=19/2,w=21/2
c
c  Basilio 22 Junio 1994
c  Basilio 23 Junio 1994 modificacion oam semienteros
c ...................................................................

        subroutine leelineasii(filelineas,iln,atom,istage,wlengt,zeff,energy,
     &                       loggf,mult,design,tam,alfa,sigma)
        integer mult(2)
        real tam(2),loggf,wlengt,zeff,alfa,sigma
        character design(2)*1,atom*2,linea*100  
        character*1 code1(21)
        character*200 nomlineas
        character*200 filelineas
        character*33 mensaje
c       character*20 cadena
        common/ficlineas/nomlineas
        common/numerito/num
        data code1/'p','1','f','2','h','3','k','4','m','5','o'
     & ,'6','r','7','t','8','u','9','v','0','w'/        


        ican=1
        mensaje=' containing the atomic parameters' 

        call cabecera(ican,filelineas,mensaje,ifail)
        if(ifail.eq.1)goto 999

c contamos las lineas
        num=0
        nxx=0
        do while(nxx.ne.iln)
           read(ican,'(a)',err=997)linea 
           call busco('=',linea,1,igu) 
           read (linea(1:igu-1),*,err=998) nxx
           num=num+1
        end do
        close(ican)
        call cabecera(ican,filelineas,mensaje,ifail)

c saltamos las num-1 primeras lineas

        do i=1,num-1          
           read(ican,'(a)',err=997)linea
        end do

c leemos la linea como una cadena

        read(ican,'(a)',err=998)linea
c	print*,'leelineasii busco',iln,'encuentro',linea
        close(ican)


c        print*,'leelineasii :',linea

c buscamos el '=' y leemos nxx y atom
        call busco('=',linea,1,igu) 
        read (linea(1:igu-1),*,err=998) nxx
        atom=linea(igu+1:igu+3)

c buscamos el proximo blanco 
        call busco(' ',linea,igu,ibl1) 

        istage=nteroleido(ibl1,ibl2,linea)

c	print*,'istage=',istage
        if(istage.gt.2)then
            print*,' '
            print*,'STOP: No ionization stages higher than 2 can be handled'
            write(*,100) 'The ionization stage',istage,'has been found in the file containing the atomic parameters'
            print*,' '
            print*,'_______________________________________________________________________________'
            stop
        endif
        wlengt=realleido(ibl2,ibl3,linea)
c        print*,'wave=',wlengt,ibl2,ibl3
        zeff=realleido(ibl3,ibl30,linea)
        energy=realleido(ibl30,ibl4,linea)
        loggf=realleido(ibl4,ibl5,linea)
        
c	print*,'wave=',wlengt
c	print*,'zeff=',zeff
c	print*,'energy=',energy
c	print*,'loggf=',loggf

c buscamos el proximo '-'. 
        call busco('-',linea,ibl5,iguion) 

c buscamos si existe un corchete antes del guion.El entero esta a la izquierda del corchete
        call busco('[',linea,ibl5,iopen1) 
        if(iopen1.lt.iguion)then
           read (linea(iopen1-2:iopen1-1),'(i2)',err=998)mult(1)
           call busco('/',linea,iopen1,ibar1) 
           read (linea(iopen1+1:ibar1-1),'(i2)',err=998)ntero1
           design(1)=code1(ntero1)
        else             !no tenemos corchete a la izquierda del -
c El entero empieza 7 characteres antes del guion
           call busco('-',linea,ibl5,iguion) 
           iint1=iguion-7
           read (linea(iint1:iint1+1),'(i2)',err=998)mult(1)
           ich1=iint1+2 
           design(1)=linea(ich1:ich1)
        end if
        read (linea(ich1+1:ich1+4),'(f4.1)',err=998)tam(1)


c buscamos si existe un corchete despues del guion.
c        print*,'____________________________________'
        call busco('[',linea,iguion,iopen2) 
c        print*,'iopen2=',iopen2
c        print*,'____________________________________'

        if(iopen2.lt.100)then
           read (linea(iopen2-2:iopen2-1),'(i2)',err=998)mult(2)
           call busco('/',linea,iopen2,ibar2) 
           read (linea(iopen2+1:ibar2-1),'(i2)',err=998)ntero2
           design(2)=code1(ntero2)
        else             !no tenemos corchete a la derecha del -
c El entero empieza 7 characteres despues del guion
           iint2=iguion+1
           read (linea(iint2:iint2+1),'(i2)',err=998)mult(2)
           ich2=iint2+2 
           design(2)=linea(ich2:ich2)  
        end if
        read (linea(ich2+1:ich2+4),'(f4.1)',err=998)tam(2)

c        print*,'en leelineas2',mult(1),design(1),tam(1)
c        print*,'en leelineas2',mult(2),design(2),tam(2)
        

        alfa=realleido(ich2+5,ial1,linea)
        sigma=realleido(ial1,ial2,linea)

c        print*,'alfa,sigma=',alfa,sigma
c        print*,'leelineasii num:',num
c        print*,'leelineasii iln:',iln

        return

999     print*,' '
        print*,'STOP: The file containing the atomic parameters does NOT exist.'
        print*,' '
        print*,'_______________________________________________________________________________'
        stop

998     print*,' '
        print*,'STOP: The file containing the atomic parameters is not written correctly'
        print*,'_______________________________________________________________________________'
        stop

997     print*,' '
        write(*,110) 'STOP: The line ',iln,'appearing in the profiles and/or wavelength'
        print*,'grid does NOT exist in the file containing the atomic parameters.'
        print*,' '
        print*,'_______________________________________________________________________________'
        stop



100     format(1x,a20,i2,1x,a59)
110     format(1x,a15,i2,1x,a43)


        end

c ___________________________________________________________________
c realleido      -lee un real si esta entre dos blancos.
c
c entradas:
c     ini        - es la posicion de un blanco anterior al real 
c     linea      - es la cadena a leer
c salidas:
c     ifi        - es la posicion del primer blanco posterior al real
c     realleido  - es el real buscado 
c
c Basilio -22 de Junio 1994
c ...................................................................

        real*4 function realleido(ini,ifi,linea)
        character linea*100 
	character str*13
   
	str='             '
c buscamos el primer 'no blanco' siguiente a ini 
        call nobusco(' ',linea,ini,inbl1) 
c buscamos el proximo 'blanco'  
        call busco(' ',linea,inbl1,ifi)
c leemos el real contenido en linea(inbl1:ifi-1)        
   	str(1:ifi-inbl1)=linea(inbl1:ifi-1)
        read (str,'(f13.5)',err=998) realleido     
c        print*,'inibl1=',inbl1,ifi-1
c        read (linea(inbl1:ifi-1),'(f13.5)',err=998) realleido
c	Esto lo tengo que hacer, en lugar de la linea anterior,
c	porque el puto compilador del KIS no le da la gana de 
c	cargar linea(inbl1:ifi-1) en realleido.

c	open(1,file='gilipollas')
c	write(1,*) linea(inbl1:ifi-1)
c	close(1)
	
c	open(1,file='gilipollas')
c	read(1,*) realleido
c	close(1,status='DELETE')
		        
        return

 998    print*,' '
        print*,'STOP: The file containing the atomic parameters is not written correctly.' 
        print*,' '
        print*,'_______________________________________________________________________________'
        stop
       end
c ___________________________________________________________________
c nteroleido      -lee un entero si esta entre dos blancos.
c
c entradas:
c     ini        - es la posicion de un blanco anterior al entero 
c     linea      - es la cadena a leer
c salidas:
c     ifi        - es la posicion del primer blanco posterior al entero
c     nteroleido  - es el entero buscado 
c
c Basilio -22 de Junio 1994
c ...................................................................

        integer*4 function nteroleido(ini,ifi,linea)
        character linea*100 
	
   
c buscamos el primer 'no blanco' siguiente a ini 
        call nobusco(' ',linea,ini,inbl1) 
c buscamos el proximo 'blanco'  
        call busco(' ',linea,inbl1,ifi)
c leemos el real contenido en linea(inbl1:ifi-1)  

        read (linea(inbl1:ifi-1),'(i1)') nteroleido

        return



998     print*,' '
        print*,'STOP: The file containing the atomic parameters is not written correctly.'
        print*,' '
        print*,'_______________________________________________________________________________'
        stop
        end   

c _______________________________________________________________
c busco          -da la posicion de un character que se encuentre
c                 a partir de una posicion dada en un string  
c entradas:
c     ch         - es el character buscado  
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................
        subroutine busco(ch,linea,ini,ifi)
        character ch*1,linea*100 

        i=ini
        do while(linea(i:i).ne.ch.and.i.lt.100)
           i=i+1
        end do
        ifi=i
c       print*,ifi
        return 
        end
c _______________________________________________________________
c nobusco         -da la primera posicion en que no aparece
c                  un character dado, comenzando la busqueda
c                  a partir duna posicion dada de un string  
c entradas:
c     ch         - es el character (NO) buscado  
c     ini        - es la posicion a partir de la cual comienza la
c                  busqueda
c salidas:
c     ifi        - es la posicion buscada 
c
c Basilio -22 de Junio 1994
c ................................................................

        subroutine nobusco(ch,linea,ini,ifi)
        character ch*1,linea*100 

        i=ini
        do while(linea(i:i).eq.ch.and.i.lt.100)
           i=i+1
        end do
        ifi=i
        return 
        end




