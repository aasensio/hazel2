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
	
	integer function mreadi3(ican,iciclo,idefault)

	character espacio*80,dato*300,linea*300,dat*13,d*1
	integer n(50),idefault,ican,iciclo
	character*100 control
	common/canal/icanal
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
              mreadi3=idefault   !si no hay nada escrito pongo el default
	      goto 25
	endif
	ivar=ivar-1
	if(iciclo.ge.ivar)then
	   mreadi3=n(ivar)
	else
           mreadi3=n(iciclo)
	end if

25	print*,espacio(1:lcom),mreadi3
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

	real*4 function mreadr3(ican,iciclo,rdefault)
	integer ican,iciclo
	real*4 x(50),rdefault
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
              mreadr3=rdefault   !si no hay nada escrito pongo el default
	      goto 26
	endif

	ivar=ivar-1
	if(iciclo.ge.ivar)then
	   mreadr3=x(ivar)
	else
           mreadr3=x(iciclo)
	end if

	if(lcom.gt.30)then
	   lcom=30
           espacio(lcom:lcom)=':'
        end if
26	write(*,100)espacio(1:lcom),mreadr3
100	format(1x,a30,1x,1pd9.2)
	return
	end
c _____________________________________________________________________________
