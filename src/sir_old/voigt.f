C*****************************************************************  
C This vectorizable Voigt function is based on the paper by Hui, *  
C Armstrong and Wray, JQSRT 19, 509 (1977). Errors become        *  
C significant (at the < 1% level) around the knee between the    *  
C Doppler core and the damping wings for 0.0 < a < 0.001. The    *  
C normalization is such that the integral is equal to SQRT(PI).  *  
C If J <> 0 this function returns the dispersion function.       *  
C -------------------------------------------------------------- *  
C Authors: Jack Harvey, Aake Nordlund.                           *  
C Modified by Sami Solanki (1985).                               *  
C Modified by A.D. Wittmann (1986) to include F(a,-v) and F(0,v) *  
C -------------------------------------------------------------- *  
C Last Update: 18-APR-86.                                        *  
C*****************************************************************  
      FUNCTION VOIGT(A,VV,J)
      COMPLEX Z 
      DIMENSION XDWS(28),YDWS(28)   
      DATA A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6/   
     *122.607931777104326,214.382388694706425,181.928533092181549,  
     *93.155580458138441,30.180142196210589,5.912626209773153,  
     *.564189583562615,122.60793177387535,352.730625110963558,  
     *457.334478783897737,348.703917719495792,170.354001821091472,  
     *53.992906912940207,10.479857114260399/
      DATA XDWS/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,   
     *3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./,YDWS/   
     *9.9335991E-02,1.9475104E-01,2.8263167E-01,3.5994348E-01,  
     *4.2443639E-01,4.7476321E-01,5.1050407E-01,5.3210169E-01,  
     *5.4072434E-01,5.3807950E-01,5.0727350E-01,4.5650724E-01,  
     *3.9993989E-01,3.4677279E-01,3.0134040E-01,1.7827103E-01,  
     *1.2934799E-01,1.0213407E-01,8.4542692E-02,7.2180972E-02,  
     *6.3000202E-02,5.5905048E-02,5.0253846E-02,4.1812878E-02,  
     *3.5806101E-02,3.1311397E-02,2.7820844E-02,2.5031367E-02/ 

 
      V=ABS(VV) 
      IF(A.NE.0) GOTO 1 
      IF(J.NE.0) GOTO 3 
      VOIGT=EXP(-V*V)   
      RETURN
   3  IF(V.GT.XDWS(1)) GOTO 4   
      D=V*(1.-.66666667*V*V)
      GOTO 8
   4  IF(V.GT.XDWS(28)) GOTO 5  
      K=27  
      DO 7 I=2,27   
      IF(XDWS(I).LT.V) GOTO 7   
      K=I   
      GOTO 6
   7  CONTINUE  
   6  KK=K-1
      KKK=K+1   
      D1=V-XDWS(KK) 
      D2=V-XDWS(K)  
      D3=V-XDWS(KKK)
      D12=XDWS(KK)-XDWS(K)  
      D13=XDWS(KK)-XDWS(KKK)
      D23=XDWS(K)-XDWS(KKK) 
      D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)* 
     *D1*D2/(D13*D23)   
      GOTO 8
   5  Y=.5/V
      D=Y*(1.+Y/V)  
   8  VOIGT=5.641895836E-1*D
   9  IF(VV.LT.0.) VOIGT=-VOIGT 
      RETURN
   1  Z=CMPLX(A,-V) 
      Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/   
     *(((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)  
      IF(J.NE.0) GOTO 2 
      VOIGT=REAL(Z)
      RETURN
   2  VOIGT=.5*AIMAG(Z) 

      GOTO 9
      END 
  
C*****************************************************************  
C This vectorizable Voigt function is based on the paper by Hui, *  
C Armstrong and Wray, JQSRT 19, 509 (1977). Errors become        *  
C significant (at the < 1% level) around the knee between the    *  
C Doppler core and the damping wings for 0.0 < a < 0.001. The    *  
C normalization is such that the integral is equal to SQRT(PI).  *  
C If J <> 0 this function returns the dispersion function.       *  
C -------------------------------------------------------------- *  
C Authors: Jack Harvey, Aake Nordlund.                           *  
C Modified by Sami Solanki (1985).                               *  
C Modified by A.D. Wittmann (1986) to include F(a,-v) and F(0,v) *
C convertida en rutina por Basilio Ruiz (1993)                   *
C -------------------------------------------------------------- *  
C Last Update: 22-jun 93.                                        *  
C*****************************************************************  
      subroutine VOIGT2(A,VV,H,F)
      COMPLEX Z 
      DIMENSION XDWS(28),YDWS(28)   
      DATA A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6/   
     *122.607931777104326,214.382388694706425,181.928533092181549,  
     *93.155580458138441,30.180142196210589,5.912626209773153,  
     *.564189583562615,122.60793177387535,352.730625110963558,  
     *457.334478783897737,348.703917719495792,170.354001821091472,  
     *53.992906912940207,10.479857114260399/
      DATA XDWS/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,   
     *3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./,YDWS/   
     *9.9335991E-02,1.9475104E-01,2.8263167E-01,3.5994348E-01,  
     *4.2443639E-01,4.7476321E-01,5.1050407E-01,5.3210169E-01,  
     *5.4072434E-01,5.3807950E-01,5.0727350E-01,4.5650724E-01,  
     *3.9993989E-01,3.4677279E-01,3.0134040E-01,1.7827103E-01,  
     *1.2934799E-01,1.0213407E-01,8.4542692E-02,7.2180972E-02,  
     *6.3000202E-02,5.5905048E-02,5.0253846E-02,4.1812878E-02,  
     *3.5806101E-02,3.1311397E-02,2.7820844E-02,2.5031367E-02/


	ivsigno=1
        if(vv.lt.0)ivsigno=-1
        v=ivsigno*vv

        IF(A.eq.0)then 
	  v2=V*V       				
          H=EXP(-v2)   
      				
          IF(V.GT.XDWS(1)) GOTO 4   
          D=V*(1.-.66666667*v2)
          GOTO 8
   4      IF(V.GT.XDWS(28)) GOTO 5  
          K=27  
          DO 7 I=2,27   
             IF(XDWS(I).LT.V) GOTO 7   
             K=I   
             GOTO 6
   7      CONTINUE  
   6      KK=K-1
          KKK=K+1   
          D1=V-XDWS(KK) 
          D2=V-XDWS(K)  
          D3=V-XDWS(KKK)
          D12=XDWS(KK)-XDWS(K)  
          D13=XDWS(KK)-XDWS(KKK)
          D23=XDWS(K)-XDWS(KKK) 
          D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)* 
     *    D1*D2/(D13*D23)   
          GOTO 8
   5      Y=.5/V
          D=Y*(1.+Y/V)  
   8      F=ivsigno*5.641895836E-1*D

c si el damping no es nulo
	else

           Z=CMPLX(A,-V) 
           Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/   
     *       (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)  
           H=REAL(Z) 
           F=.5*ivsigno*AIMAG(Z) 
             
	end if

	return
        END   
