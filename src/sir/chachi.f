
	subroutine chachi(num,cha,nc)

c	nc: Dimension de "cha" en la rutina que llama a "chachi".
c	"cha" se rellena por la izquierda con blancos.
c	La rutina chachi convierte el entero num en un character cha	

	character*(*) cha
	integer num,k,j,nc

	if (nc.gt.10) then   !Esto no sirve para nada
	   print*,'El numero de caracteres debe ser menor que 10'
	   stop
	end if
        
	
        kk=10**nc 
	k=num-(num/kk)*kk      !k es el entero formado por las ultimas
			       !nc cifras de num
	do i=nc-1,0,-1
	   j=10**i
	   n=k/j
	   k=k-j*n
	   cha(nc-i:nc-i)=char(n+48)
	end do

	return
	end
c *****************************************************************************

	subroutine quitaex0(cha) !quita la extension detras del punto

	character*(*) cha
	character*100 chacho

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end
c *****************************************************************************

	subroutine quitaex(cha)

	character*(*) cha
	character*100 chacho

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end



*******************************************************************************
	subroutine quitaex2(cha,ixt)

	integer ixt
	character*(*) cha
	character*100 chacho
        

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.'.and.cha(i:i).ne.'_')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do

	if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
	else
   	     ixt=0
	endif
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end

*******************************************************************************************************************************************************
	subroutine quitaex3(cha,ixt)

	integer ixt
	character*(*) cha
	character*100 chacho
        

	long=len(cha)

	i=1

	do while(i.lt.long.and.cha(i:i).ne.'.')
	   chacho(i:i)=cha(i:i)
	   i=i+1
	end do

	if(cha(i:i).eq.'_')then
             ixt=ichar(cha(i+1:i+1))-48
	else
   	     ixt=0
	endif
	
	do j=i,long
	   chacho(j:j)=' '
	end do

	cha(1:long)=chacho(1:long)

	return
	end
*******
		   	   
	subroutine concatena(cha1,cha2,cha)

	character*(*) cha1
	character*(*) cha2
	character*(*) cha

	long1 = len(cha1)
	long2 = len(cha2)
	long3 = len(cha)

	i=1
	do while(i.lt.long1.and.cha1(i:i).ne.' ')
	   i=i+1
	end do
	if(i .eq. long1 .and. cha1(i:i).ne.' ')i=i+1
	n1=i-1

	if (long3.le.n1) then !tampoco hace nada (dimens. se definen en el principal)
	   print*,'No puedo concatenar dos nombres en uno mas corto '
	   print*,'que el primero'
	   return
	end if

	
	i=1
	do while(i.lt.long2.and.cha2(i:i).ne.' ')
	   i=i+1
	end do
	if(i .eq. long2 .and. cha2(i:i).ne.' ')i=i+1
	n2=i-1
       
	cha(1:n1)=cha1(1:n1)
	if (n1+n2.le.long3) then
	   cha(n1+1:n1+n2) = cha2(1:n2)
	   cha(n1+n2+1:long3) = ' ' 
	else
	   cha(n1+1:long3) = cha2(1:long3-n1)
	   print*,'Hemos cortado la segunda cadena de caracteres: ',cha
	end if

	return
	end
