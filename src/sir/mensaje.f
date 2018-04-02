
	subroutine mensaje(numero,mensaje1,mensaje2,mensaje3)

	character*80 mensaje1,mensaje2,mensaje3
	integer numero,icanal
	common/canal/icanal

        print*,' '
        write(*,*) mensaje1
	if (numero.ge.2) write(*,*) mensaje2
	if (numero.ge.3) write(*,*) mensaje3
        print*,' '
        print*,'_______________________________________________________________________________'


c	open(icanal,file=control,fileopt='eof')
c        write(icanal,*) ' '
c        write(icanal,*) mensaje1
c	if (numero.ge.2) write(icanal,*) mensaje2
c	if (numero.ge.3) write(icanal,*) mensaje3
c        write(icanal,*) ' '
c        write(icanal,*) '_______________________________________________________________________________'
c	close(icanal)
        stop 

	end
