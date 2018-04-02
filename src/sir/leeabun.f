c El fichero de abundancias debe contener una cabecera
c el numero de lineas de la cabecera no debe exceder de 20 y
c cada linea no debe tener mas de 70 caracters ( si tiene
c mas no se escribiran al modificar e fichero)
c La cabecera debe acabar necesariamente con una linea que utilice alguno
c de los simbolos (ver el 'if' de lectura en el 'do while')

	subroutine leeabun(iesc,abu)           !iesc=1 escribo

        parameter (na=92)
	character mensaje*26
        integer in
        real aa,abu(na)   
        character*100 fichabun
        common/fichabun/fichabun
	common/calerr/icalerr     !si calculo errores=1 else =0 

	mensaje=' containing the abundances' 
	ican=1
	call cabecera(ican,fichabun,mensaje,ifail)
	if(ifail.eq.1)goto 999
	
         do i=1,na
           read(ican,*,err=100)in,aa
           abu(in)=aa          !10.**(aa-12.)
         end do 
100	 close(1)
         return


999	print*,' '
	print*,'STOP: The file containing the abundances does NOT exist:'
	print*,fichabun
	print*,' '
	print*,'________________________________________________________'
	stop

        end

