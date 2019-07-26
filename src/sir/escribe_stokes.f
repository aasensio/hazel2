        subroutine escribe_stokes(stok)
        include 'PARAMETER'  !por kl, kld
         
        real*4 stok(*)
c        real*4 si(kld),sq(kld),su(kld),sv(kld)
        integer ist(4),i,k,ntot
	integer ntl,nlin(kl),npas(kl),nble(kl)
	real*4 dlamda(kld)
	character*100 Stokesfilename
        common/Malla/ntl,nlin,npas,nble,dlamda  !common para StokesFRsub
        common/OutputStokes/Stokesfilename
                  
        ican=2 
         
	open(ican,file=Stokesfilename)
c contamos el numero de puntos	
	ntot=0
	do i=1,ntl
           do j=1,npas(i)
	      ntot=ntot+1
c	      print*,'escribe_stokes.f',stok(ntot)
	    end do
	end do


	k=0
	do i=1,ntl
           do j=1,npas(i)
	      k=k+1
              write(ican,993)nlin(i),dlamda(k),stok(k),stok(k+ntot),stok(k+2*ntot),stok(k+3*ntot)
           end do
        end do
        
        close(2)
        return
        
c formato de escritura
993     format(1x,i5,1x,f11.4,1x,4(e14.7,1x))       
        
        end
c______________________________________________________________________

        