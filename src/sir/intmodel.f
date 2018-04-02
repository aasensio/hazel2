c esta rutina interpola en tau el modelo ttau,tt,ppe...y da la salida en
c ttau,tt,ppe...por medio de un polinomio de grado ngrado

	subroutine intmodel(ngrado,n,tau,ttau,tt,ppe,mmic,hh,vvz,gg,ffi)

	include 'PARAMETER'  !por kt
	real ttau(*),tt(*),ppe(*),hh(*),MMIC(*),VVZ(*),gg(*),ffi(*)
	real tau(*),t(kt),pe(kt),h(kt),vz(kt),MIC(kt),G(kt),FI(kt)
	real xa(11),ya(11)

c interpolaremos las presiones en logaritmos neperianos
	num=n
	do i=1,num
	   ppe(i)=alog(ppe(i))
	end do

c interpolamos
	n2=int(ngrado/2)
	
	do i=1,n
	       CALL LOCATE(TTAU,NUM,TAU(I),J)
	       n3=j-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	       do k=1,ngrado+1
		     xa(k)=ttau(n3+k)
	       end do
	       do k=1,ngrado+1
		     ya(k)=tt(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),T(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=ppe(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),pe(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=hh(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),h(I),ERROR)
	       do k=1,ngrado+1
	             ya(k)=vvz(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),vz(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=mmic(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),mic(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=gg(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),g(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=ffi(n3+k)
	       end do
               do k=2,ngrado+1
                   dya=ya(k)-ya(k-1)
                   if(abs(dya).ge.350)ya(k)=ya(k)-360*int(dya/350.)
               end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),fi(I),ERROR)
	    end do
            do i=2,n
               dfi=fi(i)-fi(i-1)
               if(abs(dfi).ge.350)fi(i)=fi(i)-360.*int(dfi/360.)
            end do 
        
	    do i=1,n
	        ttau(i)=tau(i)
	        tt(i)=t(i)
		ppe(i)=exp(pe(i))
	        mmic(i)=mic(i)
	        hh(i)=h(i)
	        vvz(i)=vz(i)
	        gg(i)=g(i)
	        ffi(i)=fi(i)
	    end do

	return
	end
