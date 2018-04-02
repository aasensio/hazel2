c meves: funcion entera que devuelve 0 si el character de entrada es invisble
c y 1 si es visible
	integer function meves(ch,n)
	character ch*(*)
	integer n

        meves=1
	i=1
	ich=ichar(ch(1:1))
	do while((ich.eq.0.or.ich.eq.32.or.ich.eq.9).and.i.lt.n)
           i=i+1 
	   ich=ichar(ch(i:i))
c           print*,i,ich
        end do 
c        print*,i
        if(i.eq.n)meves=0
        return
	end 

