c lee una linea cuyos primeros characteres son un string hasta el simbolo :
c luego vienen una
c serie de n enteros (con n menor o igual a 15).
c La separacion entre los enteros se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna el entero
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna el ultimo entero.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993
	
	integer function mreadi2(ican,iciclo)

	character espacio*80,dato*300,linea*300,dat*13,d*1
	integer n(50),ican,iciclo
	character*100 control
c	common/canal/icanal
	common/nombrecontrol/control

	read (ican,'(a)') linea
	long=len(linea)

        i=1
	do while(linea(i:i).ne.':'.and.i.lt.long)
           i=i+1
        end do
        lcom=i 

	i=lcom+1
	do while(linea(i:i).ne.'!'.and.i.lt.long)
           i=i+1
        end do
        long=i 
  
	espacio(1:lcom)=linea(1:lcom)
	leng=long-lcom
	dato(1:leng-1)=linea(lcom+1:long-1)
	dato(leng:leng)=';'
	ivar=1
	ii=0 
	dat(1:13)='             '  
            
	do 11 i=1,leng
	   d=dato(i:i)
	   if(d.eq.'.')then
             print*,' '
	     print*,'STOP: An integer is expected in one of the fields of the control file.'
             print*,'      A real has been found instead! '
             print*,'_______________________________________________________________________________'
c	     open(icanal,file=control,fileopt='eof')	
c	     write(icanal,*) '_______________________________________________________________________________'
c             write(icanal,*) ' '
c	     write(icanal,*) 'STOP: An integer is expected in one of the fields of the control file.'
c             write(icanal,*) '      A real has been found instead! '
c	     write(icanal,*) '_______________________________________________________________________________'
c	     close(icanal)
	     stop
	   end if	     
           if(d.eq.' ')goto 11
           if(d.ne.','.and.d.ne.';')then
              ii=ii+1
	      dat(ii:ii)=d
	   else
              if(ii.eq.0)goto 11                
              ii=0
              if(dat(1:1).eq.'*')then
                 n(ivar)=1000
              else
                read (dat,'(i13)') n(ivar)
              end if
	      ivar=ivar+1
	      dat(1:13)='             '      
	   end if
11	continue
	if(ivar.eq.1)then 
              mreadi2=0   !si no hay nada escrito lo hago cero. LRB
	      goto 25
	endif
	ivar=ivar-1
	if(iciclo.ge.ivar)then
	   mreadi2=n(ivar)
	else
           mreadi2=n(iciclo)
	end if

25	print*,espacio(1:lcom),mreadi2
	return
	end
c _____________________________________________________________________________
c lee una linea cuyos primeros characteres son un string (hasta el simbolo : )
c y luego vienen una
c serie de n real*4 (con n menor o igual a 15).
c La separacion entre los real*4 se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna el real*4
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna el ultimo real*4.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993

	real*4 function mreadr2(ican,iciclo)
	integer ican,iciclo
	real*4 x(50)
	character espacio*80,dato*300,linea*300,dat*13,d*1

	read (ican,'(a)') linea
	long=len(linea)

        i=1
	do while(linea(i:i).ne.':'.and.i.lt.long)
           i=i+1
        end do
        lcom=i 

        i=lcom+1
	do while(linea(i:i).ne.'!'.and.i.lt.long)
           i=i+1
        end do
        long=i 
  
	espacio(1:lcom)=linea(1:lcom)
	leng=long-lcom
	dato(1:leng-1)=linea(lcom+1:long-1)
	dato(leng:leng)=';'
	ivar=1
	ii=0 
	dat(1:13)='             '      
	do 11 i=1,leng
	   d=dato(i:i)
           if(d.eq.' ')goto 11 
           if(d.ne.','.and.d.ne.';')then
              ii=ii+1
	      dat(ii:ii)=d
	   else
              if(ii.eq.0)goto 11
c si es un entero lo pasamos a real
              j=1 
              do while(dat(j:j).ne.'.'.and.dat(j:j).ne.'e'.and.
     &                 dat(j:j).ne.'d'.and.j.le.ii)
                 j=j+1
              end do 
              if(j.gt.ii)dat(ii+1:ii+1)='.'
                                 
              ii=0
              read (dat,'(f13.5)') x(ivar)
	      ivar=ivar+1
	      dat(1:13)='             '      
	   end if
11	continue

	if(ivar.eq.1)then 
              mreadr2=0.   !si no hay nada escrito lo hago cero. LRB
	      goto 26
	endif

	ivar=ivar-1
	if(iciclo.ge.ivar)then
	   mreadr2=x(ivar)
	else
           mreadr2=x(iciclo)
	end if

	if(lcom.gt.30)then
	   lcom=30
           espacio(lcom:lcom)=':'
        end if
26	write(*,100)espacio(1:lcom),mreadr2
100	format(1x,a30,1x,1pd9.2)
	return
	end
c _____________________________________________________________________________
c lee una linea cuyos primeros characteres son un string (hasta el simbolo :) 
c y luego vienen una
c serie de n cadenas o palabras (con n menor o igual a 15).
c La separacion entre las cadenas se supone realizada por los simolos , o ;
c Se ignoran los blancos 
c Esta function retorna la cadena
c que ocupa el orden iciclo (el icilcloesimo vamos), en caso de que iciclo sea
c mayor que n retorna la ultima cadena.
c Ignora todo aquello que se encuentre a la derecha del simbolo !
c Basilio -Enero 1993
	character*100 function mreadc2(ican,iciclo)
	integer ican,iciclo
	character*100 x(100)
	character titulo*100,dato*400,linea*400,dat*100,d*1,titprint*150

	read (ican,'(a)') linea
	long=len(linea)

        i=1
	do while(linea(i:i).ne.':'.and.i.lt.long)
           i=i+1
        end do
        lcom=i 

        i=lcom+1
	do while(linea(i:i).ne.'!'.and.i.lt.long)
           i=i+1
        end do
        long=i 
  
	titulo(1:lcom)=linea(1:lcom)
	leng=long-lcom
	dato(1:leng-1)=linea(lcom+1:long-1)
	dato(leng:leng)=';'
	ivar=1
	ii=0 
	dat(1:100)='                                                                                                    '  
c	print*,'dato: ',dato(1:leng),'/'    
	do 11 i=1,leng
	   d=dato(i:i)
           if(d.eq.' ')goto 11 
           if(d.ne.','.and.d.ne.';')then
c              print*,'if1: ',d,':'       
	      ii=ii+1
	      dat(ii:ii)=d
	   else
              if(ii.eq.0)goto 11                
              x(ivar)(1:100)=dat(1:100)
              ii=0
	      ivar=ivar+1
	      dat(1:100)='                                                                                                    '  
	   end if
11	continue
	ivar=ivar-1
	if(ivar.ne.0)then 
	  if(iciclo.ge.ivar)then
	     mreadc2(1:100)=x(ivar)(1:100)
	  else
             mreadc2(1:100)=x(iciclo)(1:100)
	  end if
	else
             mreadc2(1:100)='                                                                                                    '  

        endif

	titprint=titulo(1:lcom)//' '//mreadc2(1:48)
	print*,titprint
c	print*,titulo(1:lcom),' ',mreadc2(1:48)
        
	return
	end




c _____________________________________________________________________________

	subroutine mreadmalla(ican,numblends,n,ini,paso,fin,ierror,blanco)

	include 'PARAMETER'
	real*4 x(50)
	character(len=80) dato,linea
	character(len=13) dat
	character(len=1) d
	integer n(200),blanco,ican,numblends,ierror, caca
	real*4 ini,paso,fin

	common/canal/icanal

	ierror=0
	blanco=0
	read (ican,'(a)',end=999) linea
	longg=len(linea)

	if(longg.eq.0.or.longg.eq.1)then   !la linea esta vacia, o no esta bien escrita
           blanco=1
           return
	endif

	
        
	lcom=0 

        i=lcom+1
	do while(linea(i:i).ne.':'.and.i.lt.longg)
           i=i+1
        end do
	if(i.eq.longg)then   !la linea esta vacia, o no esta bien escrita
           blanco=1
           return
	endif
        long=i 
  
	leng=long-lcom
	dato(1:leng-1)=linea(lcom+1:long-1)
	dato(leng:leng)=';'
	ivar=1
	ii=0 
	dat(1:13)='             ' 
     
	do 11 i=1,leng
	   d=dato(i:i)
	   if(d.eq.'.')then
 	     stop 'error en mreadmalla'
	   end if	     
           if(d.eq.' ')goto 11
           if(d.ne.','.and.d.ne.';')then
              ii=ii+1
	      dat(ii:ii)=d
	   else
              if(ii.eq.0)goto 11                
              ii=0			  
              read (dat,'(i13)') n(ivar)
	      ivar=ivar+1
	      dat(1:13)='             '      
	   end if
11	continue
	if(ivar.eq.1)then
	      print*,' ' 
              print*,'STOP: Incorrect format in the file containing the wavelength grid.'
	      print*,' ' 
              print*,'_______________________________________________________________________________'
c              write(icanal,*) 'STOP: Incorrect format in the file containing the wavelength grid. '
   	endif
	ivar=ivar-1
	numblends=ivar

c	Ahora leemos las lambdas:

	leng=longg-long
	dato(1:leng-1)=linea(long+1:longg-1)
	dato(leng:leng)=';'
	ivar=1
	ii=0 
	dat(1:13)='             ' 

	do 12 i=1,leng
	   d=dato(i:i)
           if(d.eq.' ')goto 12 
           if(d.ne.','.and.d.ne.';')then
              ii=ii+1
	      dat(ii:ii)=d
	   else
              if(ii.eq.0)goto 12
c             si es un entero lo pasamos a real:
              j=1 
              do while(dat(j:j).ne.'.'.and.dat(j:j).ne.'e'.and.
     &                 dat(j:j).ne.'d'.and.j.le.ii)
                 j=j+1
              end do 
              if(j.gt.ii)dat(ii+1:ii+1)='.'
                                 
              ii=0
c		print*,'entro',dat
              read (dat,'(f13.5)') x(ivar)
	      ivar=ivar+1
	      dat(1:13)='             '      
	   end if
12	continue

	if(ivar.gt.4)then 
           print*,' '
           print*,'WARNING: There are more than the three items required to determine the wavelength'
           print*,'range in the file containing the wavelength grid. The first three are adopted'
	   print*,' '
           print*,'_______________________________________________________________________________'
c           write(icanal,*) 'WARNING: There are more than the three items required to determine the wavelength'
c           write(icanal,*) 'range in the file containing the wavelength grid. The first three are adopted.'
	endif
	if(ivar.lt.4)then  !falta algun numero 
           print*,' '
           print*,'STOP: Some parameter required to determine the wavelength range is missing in'
	   print*,'      the file containing the wavelength grid. Check the values after the : symbol.'
	   print*,' '
           print*,'_______________________________________________________________________________'
c           write(icanal,*) 'STOP: Some parameter required to determine the wavelength range is missing in'
c           write(icanal,*) '      the file containing the wavelength grid. Check the values after the : symbol.'
	   stop
	endif
	  

	ini=x(1)
	paso=x(2)
	fin=x(3)



	return

999	ierror=1
	return
	end
