        subroutine planck2(temp,x,bratio,bp,bt)

        real*4 expo,temp
        real*4 bp,bt,bratio

        a=1.43880/x
        c=(1.1910627e-5/x)/x**4

        bp=dplnck(temp,x) 

        if(bratio.eq.1.)then
           bt=dtplanck(temp,x) 
           return
        else
           expo=c/bp+1.
           bp=c/(bratio*expo-1.)
           bt=bp*bp*bratio*a*expo/c/temp/temp
        end if
	return     

        end   
