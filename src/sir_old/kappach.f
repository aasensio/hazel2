c--------------------------------------------------------------------
c	la subrutina kappach es una modificacion dela function kappac
c	para calcular las derivadas del coef. de absorcion (kac) con la
c	temperatura (kat), y la presion electronica (kap)
c
      	subroutine kappach(lambda,tt,ppe,p,dp,ddp,kac,kat,kap)

	real*4 kac,kat,kap
	real p(10),dp(10),ddp(10),pg(99),dpg(99),ddpg(99)
	real dz3(12)
	integer ivar(10)

      real mgatom,nair		!,wavelt/0./
      real*4 dfreq,lnfreq,dobtheta,dobdtheta,lambda
      logical wrepet
      dimension z3(12)
      real ghel(12)/1.,3.,1.,9.,9.,3.,3.,3.,1.,9.,20.,3./,
     *chihel(12)/0.,19.819,20.615,20.964,20.964,21.217,21.217,22.718,22.
     *920,23.007,23.073,23.087/,g(12),tcatom(27)/3.9999e3,5e3,6e3,7e3,8e
     *3,9e3,1e4,1.1e4,1.2e4,1.3e4,1.4e4,1.5e4,1.6e4,1.7e4,1.8e4,1.9e4,2e
     *4,2.1e4,2.2e4,2.3e4,2.4e4,2.5e4,2.6e4,2.7e4,2.8e4,2.9e4,3.00001e4/
     *,fcatom(27)/1.728,6.660,9.964,12.340,14.134,15.540,16.676,17.614,1
     *8.402,19.076,19.662,20.172,20.626,21.030,21.394,21.722,22.022,22.2
     *96,22.548,22.78,22.996,23.196,23.384,23.56,23.724,23.878,24.024/,
     *fsodiu(27)/1.638,3.914,5.461,6.586,7.448,8.130,8.686,9.152,9.546,
     *9.887,10.184,10.448,10.682,10.894,11.084,11.258,11.418,11.566,11.7
     *04,11.830,11.950,12.062,12.166,12.264,12.358,12.448,12.532/
      real fmatom(27)/1.784,5.268,7.618,9.318,10.61,11.628,12.452,1
     *3.134,13.71,14.204,14.632,15.008,15.342,15.64,15.906,16.15,16.37,1
     *6.574,16.76,16.934,17.094,17.244,17.384,17.514,17.638,17.754,17.86
     *6/,wcatom(15)/0.,1.1005e-5,1.1005e-5,1.2395e-5,1.2395e-5,1.4445e-5
     *,1.4445e-5,2.178e-5,2.912e-5,3.4625e-5,3.4625e-5,4.7295e-5,
     *4.7295e-5,5.4785e-5,7.1e-5/,wsodiu(21)/0.,1.2e-5,1.5e-5,1.6e-5,1.7
     *e-5,1.8e-5,1.95e-5,2.1e-5,2.2e-5,2.3e-5,2.4125e-5,2.4125e-5,2.75e-
     *5,4.0845e-5,4.0845e-5,6.366e-5,8.1455e-5,8.1455e-5,8.9455e-5,8.945
     *5e-5,9.85e-5/,wmatom(20)/1.3999e-5,1.5e-5,1.6215e-5,1.6215e-5,
     *2.5135e-5,2.5135e-5,2.75e-5,3.7565e-5,3.7565e-5,4.8845e-5,
     *4.8845e-5,6.5495e-5,6.5495e-5,7.2345e-5,7.2345e-5,7.2915e-5,
     *7.2915e-5,8.1135e-5,8.1135e-5,9.e-5/
      real ccatom(15,9)/
     *-2.75,-3.02,-1.16,-1.17,1.33,1.32,8.05,7.79,7.66,7.68,8.12,8.07,
     *8.28,8.44,8.14,-2.75,-3.02,-1.48,-1.49,0.66,0.65,5.99,5.73,5.61,
     *5.63,5.86,5.79,5.95,5.97,5.67,-2.75,-3.02,-1.68,-1.69,0.22,0.20,
     *4.58,4.32,4.21,4.20,4.36,4.25,4.36,4.35,4.02,-2.75,-3.02,-1.83,
     *-1.84,-0.10,-0.11,3.56,3.30,3.20,3.18,3.27,3.10,3.19,3.08,2.83,
     *-2.75,-3.02,-1.94,-1.95,-0.33,-0.34,2.78,2.52,2.38,2.36,2.43,2.22,
     *2.34,2.25,1.93,-2.75,-3.02,-2.08,-2.09,-0.65,-0.67,1.67,1.41,1.27,
     *1.20,1.24,1.03,1.08,0.97,0.66,-2.75,-3.02,-2.27,-2.29,-1.10,-1.12
     *,0.11,-0.14,-0.21,-0.42,-0.40,-0.64,-0.64,-0.78,-1.09,-2.75,-3.02
     *,-2.37,-2.38,-1.38,-1.41,-0.75,-0.98,-1.16,-1.29,-1.28,-1.57,-1.55
     *,-1.70,-2.01, -2.75,-3.02,-2.51,-2.53,-1.88,-1.91,-1.70,-1.91,-2.1
     *0,-2.24,-2.24,-2.55,-2.54,-2.70,-3.01/ ! -(log(kappa)+20)
      real cmatom(20,8)/
     *-0.57,-1.55,-2.08,-0.20,-0.75,2.14,2.01,2.02,3.03,2.74,2.79,2.49,2
     *.71,2.59,2.75,2.74,3.75,3.67,3.78,3.61,-0.87,-1.61,-2.10,-0.89,-1.
     *44,1.01,0.88,0.88,1.56,1.26,1.30,1.00,1.19,1.08,1.23,1.22,2.05,1.9
     *4,1.98,1.87,-1.20,-1.71,-2.14,-1.34,-1.89,0.25,0.14,0.07,0.57,0.27
     *,0.32,0.00,0.17,0.07,0.20,0.18,0.88,0.78,0.82,0.69,-1.72,-1.99,
     *-2.30,-1.92,-2.47,-0.74,-0.85,-0.98,-0.69,-0.97,-0.95,-1.22,-1.12,
     *-1.23,-1.10,-1.12,-0.59,-0.70,-0.68,-0.81,-2.05,-2.24,-2.48,-2.26,
     *-2.82,-1.36,-1.47,-1.65,-1.45,-1.68,-1.68,-2.04,-1.90,-2.02,-1.91,
     *-1.92,-1.49,-1.60,-1.60,-1.73,-2.53,-2.67,-2.83,-2.74,-3.30,-2.24,
     *-2.35,-2.59,-2.48,-2.73,-2.73,-3.10,-3.00,-3.12,-3.02,-3.04,-2.77,
     *-2.88,-2.88,-3.00,-2.82,-2.92,-3.06,-3.01,-3.55,-2.73,-2.82,-3.11,
     *-3.03,-3.34,-3.33,-3.66,-3.58,-3.69,-3.63,-3.64,-3.43,-3.56,-3.56,
     *-3.69,-3.18,-3.26,-3.36,-3.33,-3.85,-3.26,-3.36,-3.67,-3.63,-3.94,
     *-3.93,-4.28,-4.22,-4.34,-4.29,-4.30,-4.17,-4.30,-4.30,-4.43/ ! -log-20
      real csodiu(21,8)/8.70,  
     *8.78,8.99,9.13,9.41,10.11,10.62,9.95,9.31,9.00,8.87,10.22,9.96,9.2
     *6,11.45,10.76,10.43,10.95,10.82,11.32,11.18,8.70,8.78,8.99,9.13,9.
     *37,9.93,10.10,9.72,9.23,8.96,8.82,9.67,9.43,8.73,10.52,9.83,9.48, 
     *9.92,9.80,10.18,10.04,8.70,8.78,8.99,9.12,9.31,9.72,9.72,9.47,9.11
     *,8.88,8.75,9.31,9.07,8.34,9.86,9.18,8.84,9.24,9.12,9.41,9.27,8.70,
     *8.75,8.94,9.01,9.15,9.33,9.25,9.06,8.85,8.68,8.57,8.84,8.6,7.93,
     *9.03,8.37,8.02,8.33,8.20,8.41,8.28,8.65,8.72,8.86,8.91,9.00,9.06,
     *8.95,8.78,8.63,8.49,8.39,8.55,8.31,7.63,8.52,7.86,7.52,7.77,7.64,7
     *.79,7.66,8.60,8.61,8.62,8.63,8.62,8.59,8.47,8.33,8.23,8.13,8.02,8.
     *10,7.87,7.19,7.79,7.12,6.87,6.94,6.82,6.93,6.78,8.50,8.48,8.42,8.3
     *8,8.34,8.29,8.17,8.06,7.96,7.87,7.79,7.83,7.61,6.98,7.35,6.70,6.36
     *,6.49,6.36,6.43,6.29,8.40,8.26,8.10,8.04,7.98,7.89,7.78,7.67,7.59,
     *7.52,7.44,7.45,7.24,6.62,6.88,6.23,5.88,5.96,5.83,5.87,5.73/ ! -log-10
c     dimension acl(4)/1.38917648e+00,7.65012349e-01,-2.78423595e-01,
c    *4.04414569e-02/,bcl(4)/-2.96225925e+01,-1.45062313e+00,1.72112600e
c    *+00,-3.01692989e-01/,ccl(4)/1.19071498e+02,3.11490334e+01,
c    *-1.59590407e+01,2.34426531e+00/,dcl(4)/-1.77615290e+02,-5.63832483
c    *e+01,2.48188080e+01,-3.35147209e+00/

	ivar(1)=1
	ivar(2)=2
	ivar(3)=6
	ivar(4)=11
	ivar(5)=12
	ivar(6)=86
        ivar(7)=89
	ivar(8)=90
	ivar(9)=91
	ivar(10)=93

	do i=1,99
	   pg(i)=0.0
	   dpg(i)=0.0
	   ddpg(i)=0.0
	end do

	do i=1,10
	   j=ivar(i)
	   pg(j)=p(i)
	   dpg(j)=dp(i)
	   ddpg(j)=ddp(i)
	end do
	t=tt
	pe=ppe
	wavelt=0.	!lo introduzco nuevo
      theta=5040./t
	dtheta=-5040./(t*t)
      dobtheta=dble(theta)
	dobdtheta=dble(dtheta)
      temp25=t**2.5
      wrepet=.false.
      if(abs(lambda-wavelt).lt.1.e-8) wrepet=.true.
      if(wrepet) goto 3
      x10000=1.e+4*lambda   ! microns
      nair=1.0004
      if(x10000.lt.0.18) goto 23
      nair=refrax(x10000,15.,760.,0.)
23    freq=2.997925e+10/(lambda*nair)
      dfreq=dble(freq)
      lnfreq=alog(dfreq)
      z=alog10(freq)
      x10001=lambda**3
      x10002=1.e+8*lambda   ! angstroms
      deltak=911.3/(x10002*nair)
      ephot=deltak
      ey=ephot**(0.43+0.6*alog10(ephot+10.))
      divi1=1.+(5.9856e-2-(3.4916e-4*ey)/(1.+1.e-2*ey))*ephot**.83333333
      m0=int(sqrt(1./ephot))+1
      l0=max0(m0,4)
      t1=1.e-8/lambda
      t2=t1**2
      t3=t2**2
      x10003=1.4388/lambda
      x10004=1.e+2*lambda
      x10005=1.e+5*lambda
      x10006=1.e-34*x10002**2
      scat1=t3*(5.799e-13+1.422e-6*t2+2.784*t3)
      scat2=t3*(8.14e-13+1.28e-6*t2+1.61*t3)
      scat3=5.484e-14*t3*(1.+2.44e5*t2+5.94e-10*t2/(x10002**2-2.9e5))**2
c_______________________________________________________________________
3     hminus=0. ! *** h- ***
        dhminus=0.
	ddhminus=0.
      f1=x10003/t
	df1=-f1/t
      z1=1.-exp(-f1)
	dz1=df1*exp(-f1)
      	if(f1.lt.1.e-3)then
 	z1=f1
	dz1=df1
	end if

      if(lambda.gt.1.64189e-4) goto 5
      if(wrepet) goto 37
      if(lambda.gt.1.42e-4) goto 69
      x=x10005
      cbfree=1.e-17*(6.80133e-3+x*(1.78708e-1+x*(1.6479e-1-x*(2.04842e-2
     *-5.95244e-4*x))))
      goto 37
69    	x=16.419-x10005
      cbfree=1.00000e-17*x*(2.69818e-1+x*(2.2019e-1-x*(4.11288e-2-2.7323
     *6e-3*x)))
37    	hminus=cbfree*pg(93)*z1/pe
	dhminus=hminus*(dpg(93)+dz1/z1)
	ddhminus=hminus*(ddpg(93)-1./pe)

c	.	.	.	.	.	.	.	.	
c5     	hminus=(hminus+1.e-26*(5.3666e-3+theta*(2.7039e-2*theta-1.1493e-2)
c     *+x10004*(-3.2062+theta*(11.924-5.939*theta)+x10005*(-4.0192e-1+the
c     *ta*(7.0355-theta*3.4592e-1)))))*pe*pg(1)
c	.	.	.	.	.	.	.	.	

5	b1=11.924-5.939*theta
	db1=-5.939*dtheta
	c1=7.0355-theta*3.4592e-1
	dc1=-3.4592e-1*dtheta
	c2=x10005*(-4.0192e-1+theta*c1)
	dc2=x10005*(dtheta*c1+theta*dc1)
	b2=x10004*(-3.2062+theta*b1+c2)
	db2=x10004*(dtheta*b1+theta*db1+dc2)
        b3=2.7039e-2*theta-1.1493e-2+b2
	db3=2.7039e-2*dtheta+db2
	b4=5.3666e-3+theta*b3
	db4=dtheta*b3+theta*db3
	hminus=(hminus+1.e-26*b4)*pe*pg(1)
	dhminus=(dhminus+1.e-26*db4)*pe*pg(1)+hminus*dpg(1)
	

c	dhminus=(dhminus+1.e-26*(dtheta*(2.7039e-2*2.*theta-1.1493e-2)+
c     *	x10004*dtheta*(11.924-5.939*2.*theta)+x10005*(dtheta*(7.0355-2.*
c     *	theta*3.4592e-1))))*pe*pg(1)+hminus*dpg(1)

	ddhminus=ddhminus*pe*pg(1)+hminus/pe+hminus*ddpg(1)
c---------------------------------------------------------------------
   4  helmin=0. ! *** he- ***

	dhelmin=0.
	ddhelmin=0.
      if(deltak.gt.0.3.or.t.lt.1.5e3.or.t.gt.1.68e4) goto 73
      if(t.lt.9.2e3) goto 71
      	a1=2.46e-4+lambda*(-1.26e+1+5.67e+6*lambda)
	b1=theta*(-5.92e-4+lambda*(5.12e+1+1.3e+7*lambda))
	c1=(theta*theta)*(1.45e-2-lambda*(8.3e+1+1.8e+6*lambda))
	da1=0.
	db1=b1*dtheta/theta
	dc1=2.*c1*dtheta/theta
	helmin=(a1+b1+c1)*1.e-26
	dhelmin=(da1+db1+dc1)*1.e-26

c      helmin=(2.46e-4+lambda*(-1.26e+1+5.67e+6*lambda)+theta*(-5.92e-4
c     *+lambda*(5.12e+1+1.3e+7*lambda)+theta*(1.45e-2-lambda*(8.3e+1+
c     *1.8e+6*lambda))))*1.e-26
      goto 72
c  71  helmin=x10006*(.49245+t*(-1.0754e-4+t*(9.5114e-9-t*2.3544e-13)))
71	a1=9.5114e-9-t*2.3544e-13
	a2=-1.0754e-4+t*a1
	a3=.49245+t*a2
	helmin=x10006*a3
	da1=-2.3544e-13
	da2=a1+t*da1
	da3=a2+t*da2
	dhelmin=x10006*da3

72      helmin=helmin*pe*pg(2)
      	dhelmin=dhelmin*pe*pg(2)+helmin*dpg(2)
	ddhelmin=ddhelmin*pe*pg(2)+helmin*ddpg(2)+helmin/pe
c_______________________________________________________________________
73    clmin=0. ! *** cl- ***

        dclmin=0.
	ddclmin=0.
c     if(deltak.lt.0.019.or.deltak.gt.0.35) go to 6
c     clmin=1.e-26*10.**(acl(1)+theta*(acl(2)+theta*(acl(3)+theta*acl(4)
c    *))+deltak*(bcl(1)+theta*(bcl(2)+theta*(bcl(3)+theta*bcl(4)))
c    *+deltak*(ccl(1)+theta*(ccl(2)+theta*(ccl(3)+theta*ccl(4)))
c    *+deltak*(dcl(1)+theta*(dcl(2)+theta*(dcl(3)+theta*dcl(4)))))))
c     clmin=clmin*pe*pg(17)
c-----------------------------------------------------------------------
   6  hneutr=0. ! *** h ***
	dhneutr=0.
	ddhneutr=0.

      if(m0.gt.12) goto 10
      x2=um(12,t)
	dx2=dum(x2,t)
      x3=um(1,t)
	dx3=dum(x3,t)
      sum=0.
	dsum=0.
      do 9 i=m0,12
	ii=i
	uuu=um(ii,t)
      z3(i)=exp(uuu-x3)/float(i**3)
	um1=dum(uuu,t)
	dz3(i)=z3(i)*(um1-dx3)
      if(wrepet) goto 91
      divi2=1.+(1.72826e-1-3.45652e-1/(ephot*i**2))*ephot**.333333333
91      g(i)=divi2/divi1   ! bound-free gaunt factors
    	sum=sum+g(i)*z3(i)
9	dsum=dsum+g(i)*dz3(i)
      z2=2.
	dz2=0.
      	if(t.gt.1.3e4)then
	 z2=1.51+3.8e-5*t
	dz2=3.8e-5
	end if
      	if(t.gt.1.62e4)then
 	z2=11.41+t*(-1.1428e-3+t*3.52e-8)
	dz2=-1.1428e-3+2.*t*3.52e-8
	end if
	a1=z2
	da1=dz2
      z2=2.08966e-2*x10001*z1/a1
c	dz2=z2*(dz1/z1-da1/a1)
	dzz2=dz1/z1-da1/a1
	dz2=z2*dzz2
	gg=gff(t,x10002)
      z4=0.5/x3*(exp(x2-x3)+exp(-x3)*(gg-1.))
	ggg=dgff(t,x10002)
	dz4=-z4*dx3/x3+.5/x3*(exp(x2-x3)*(dx2-dx3)+exp(-x3)*(-dx3*(gg-1.
     *	)+ggg))
      hneutr=z2*(sum+z4)*pg(1)
	if(abs(sum+z4).gt.1.e-30)then
	    dhneutr=hneutr*(dzz2+dpg(1)+(dsum+dz4)/(sum+z4))
	else
	    dhneutr=hneutr*(dzz2+dpg(1))
	end if
	ddhneutr=hneutr*ddpg(1)
c-----------------------------------------------------------------------
10    h2min=0. ! *** h2- ***

      if(deltak.gt.0.3.or.t.lt.1.5e3.or.t.gt.1.68e4) goto 20
      h2min=x10006*(.88967+t*(-1.4274e-4+t*(1.0879e-8-t*2.5658e-13)))*
     *pe*pg(89)
	dh2min=x10006*(-1.4274e-4+2.*t*(1.0879e-8-3.*t*2.5656e-13))*pe*
     *	pg(89)+h2min*dpg(89)
	ddh2min=h2min*(1./pe+ddpg(89))
c-----------------------------------------------------------------------
20    h2plus=0. ! *** h2+ ***

      if(lambda.lt.3.8e-5.or.lambda.gt.3e-4) goto 70
      h2plus=sngl(dexp(2.30258509d0*dobtheta*(7.342d-3-(-2.409d-15+(1.028d
     *-30+(-4.23d-46+(1.224d-61-1.351d-77*dfreq)*dfreq)*dfreq)*dfreq)*df
     *req)-3.0233d3+(3.7797d2+(-1.82496d1+(3.9207d-1-3.1672d-3*lnfreq)*l
     *nfreq)*lnfreq)*lnfreq)*1.d16)*z1*(pg(1)*pe)*(pg(86)/pg(90))/(1.380
     *54*t)  ! factor 1.d16 from boltzmann constant

     	d1=sngl(2.30258509d0*dobdtheta*(7.342d-3-(-2.409d-15+(1.028d
     *-30+(-4.23d-46+(1.224d-61-1.351d-77*dfreq)*dfreq)*dfreq)*dfreq)*df
     *req))
	dh2plus=h2plus*(d1+dz1/z1+dpg(1)+(dpg(86)-dpg(90))-1./t)
	ddh2plus=h2plus*(ddpg(1)+(ddpg(86)-ddpg(90))+1./pe)
c-----------------------------------------------------------------------
70    heneut=0. ! *** he ***
	dheneut=0.
	ddheneut=0.

      if(lambda.gt.8.2610e-5) goto 30
      if(wrepet) goto 38
      if(lambda.gt.3.6788e-5) goto 36
      i0=1 ! 1 1s
      if(lambda.gt.5.0420e-6) i0=2 ! 2 3s
      if(lambda.gt.2.6003e-5) i0=3 ! 2 1s
      if(lambda.gt.3.1210e-5) i0=4 ! 2 3p odd, 2 3p odd
      if(lambda.gt.3.4210e-5) i0=6 ! 2 1p odd, 2 1p odd
      goto 38
  36  i0=8 ! 3 3s
      if(lambda.gt.6.6322e-5) i0=9 ! 3 1s
      if(lambda.gt.7.4351e-5) i0=10 ! 3 3p
      if(lambda.gt.7.8438e-5) i0=11 ! 3 1d + 3 3d
      if(lambda.gt.8.1910e-5) i0=12 ! 3 1p odd

  38  do 21 i=i0,12
      iii=i
  	pepa=ghel(i)*10.**(fkny(iii,z)-theta*chihel(i))
	dpepa=pepa*(-dtheta*chihel(i))*alog(10.)
	dheneut=dheneut+dpepa
  21  heneut=heneut+pepa
  30  if(m0.gt.12) goto 31
      sum=0.
	dsum=0.
      do 29 i=l0,12
	dsum=dsum+g(i)*dz3(i)
  29  sum=sum+g(i)*z3(i)
        if(theta.gt.3.0.or.abs(z2).lt.1.e-25.or.abs(sum+z4).lt.1.e-25)goto31
        pepo=4.*10.**(-10.992*theta)*z2*(sum+z4)
	dpepo=pepo*(-alog(10.)*10.992*dtheta+dz2/z2+(dsum+dz4)/(sum+z4))
        heneut=heneut+pepo
	dheneut=dheneut+dpepo
  31  heneut=heneut*pg(2)
	dheneut=dheneut*pg(2)+heneut*dpg(2) !correcto
	ddheneut=ddheneut*pg(2)+heneut*ddpg(2) !correcto
c-----------------------------------------------------------------------
        kac=hminus+helmin+hneutr+h2min+heneut+h2plus
      	kat=dhminus+dhelmin+dhneutr+dh2min+dheneut+dh2plus
      	kap=ddhminus+ddhelmin+ddhneutr+ddh2min+ddheneut+ddh2plus
c-----------------------------------------------------------------------
      scatt1=0.
      scatt2=0.
      scatt3=0.
	dscatt1=0.
	dscatt2=0.
	dscatt3=0.
	ddscatt1=0.
	ddscatt2=0.
	ddscatt3=0.
      if(lambda.lt.1.2e-5) goto 35
      scatt1=scat1*pg(1) ! *** scattering ***
	dscatt1=scatt1*dpg(1)
	ddscatt1=scatt1*ddpg(1)
      scatt2=scat2*pg(89)
	dscatt2=scatt2*dpg(89)
	ddscatt2=scatt2*ddpg(89)
      scatt3=scat3*pg(2)
	dscatt3=scatt3*dpg(2)
	ddscatt3=scatt3*ddpg(2)
  35  escatt=6.653e-25*pg(90)
	descatt=escatt*dpg(90)
	ddescatt=escatt*ddpg(90)
      	kac=kac+escatt+scatt1+scatt2+scatt3
	kat=kat+descatt+dscatt1+dscatt2+dscatt3
	kap=kap+ddescatt+ddscatt1+ddscatt2+ddscatt3
c-----------------------------------------------------------------------
      carbon=0. ! *** c ***
      if(t.gt.tcatom(1).and.t.lt.tcatom(27)) goto 142
      goto 99999
 142  do 143 i=2,27
      if(t.gt.tcatom(i)) goto 143
      k=i
      goto 144
 143  continue
 144  if(k.eq.27) k=26
      ksod=k
      kk=k-1
      kksod=kk
      kkk=k+1
      kkksod=kkk
      if(lambda.lt.wcatom(15)) goto 141
      goto 50
 141  y=fint(t,fcatom(kk),fcatom(k),fcatom(kkk),tcatom(kk),tcatom(k),tca
     *tom(kkk))
      if(k.gt.5) goto 45
      y=(y-fcatom(kk))/(fcatom(k)-fcatom(kk))
      l=kk
      goto 47
  45  if(k.gt.7) goto 76
      y=(y-fcatom(5))/2.542   ! (7)-(5)
      l=5
      goto 47
76    if(k.gt.12) goto 77
      y=(y-fcatom(7))/3.496   ! (12)-(7)
      l=6
      goto 47
77    if(k.gt.17) goto 46
      y=(y-fcatom(12))/1.85   ! (17)-(12)
      l=7
      goto 47
  46  y=(y-fcatom(17))/2.002   ! (27)-(17)
      l=8
  47  do 48 i=2,15
      if(lambda.gt.wcatom(i)) goto 48
      k=i
      goto 49
  48  continue
  49  kk=k-1
      y2=(lambda-wcatom(kk))/(wcatom(k)-wcatom(kk))
      y1=ccatom(kk,l)+y2*(ccatom(k,l)-ccatom(kk,l))
      y2=ccatom(kk,l+1)+y2*(ccatom(k,l+1)-ccatom(kk,l+1))
      carbon=z1*pg(6)*10.**(y*(y1-y2)-y1-20.)
	dcarbon=carbon*(dz1/z1+dpg(6))
	ddcarbon=carbon*(ddpg(6))
c-----------------------------------------------------------------------
  50  sodium=0. ! *** na ***

      if(lambda.lt.wsodiu(21)) goto 53
      goto 61
  53  y=fint(t,fsodiu(kksod),fsodiu(ksod),fsodiu(kkksod),tcatom(kksod),
     *tcatom(ksod),tcatom(kkksod))
      if(ksod.gt.3) goto 55
      y=(y-fsodiu(kksod))/(fsodiu(ksod)-fsodiu(kksod))
      l=kksod
      goto 58
  55  if(ksod.gt.5) goto 56
      y=(y-fsodiu(3))/1.987   ! (5)-(3)
      l=3
      goto 58
  56  if(ksod.gt.7) goto 74
      y=(y-fsodiu(5))/1.238   ! (7)-(5)
      l=4
      goto 58
74    if(ksod.gt.12) goto 75
      y=(y-fsodiu(7))/1.762   ! (12)-(7)
      l=5
      goto 58
75    if(ksod.gt.17) goto 57
      y=(y-fsodiu(12))/0.97   ! (17)-(12)
      l=6
      goto 58
  57  y=(y-fsodiu(17))/1.114   ! (27)-(17)
      l=7
58    do 59 i=2,21
      if(lambda.gt.wsodiu(i)) goto 59
      k=i
      goto 60
  59  continue
  60  kk=k-1
      y2=(lambda-wsodiu(kk))/(wsodiu(k)-wsodiu(kk))
      y1=csodiu(kk,l)+y2*(csodiu(k,l)-csodiu(kk,l))
      y2=csodiu(kk,l+1)+y2*(csodiu(k,l+1)-csodiu(kk,l+1))
      sodium=z1*pg(11)*10.**(y*(y1-y2)-y1-10.)
	dsodium=sodium*(dz1/z1+dpg(11))
	ddsodium=sodium*ddpg(11)
c-----------------------------------------------------------------------

  61  mgatom=0. ! *** mg ***
      if(lambda.gt.wmatom(1).and.lambda.lt.wmatom(20)) goto 62
      goto 52
  62  y=fint(t,fmatom(kksod),fmatom(ksod),fmatom(kkksod),tcatom(kksod),
     *tcatom(ksod),tcatom(kkksod))
      if(ksod.gt.3) goto 63
      y=(y-fmatom(kksod))/(fmatom(ksod)-fmatom(kksod))
      l=kksod
      goto 66
  63  if(ksod.gt.5) goto 64
      y=(y-fmatom(3))/2.992   ! (5)-(3)
      l=3
      goto 66
  64  if(ksod.gt.7) goto 78
      y=(y-fmatom(5))/1.842   ! (7)-(5)
      l=4
      goto 66
78    if(ksod.gt.12) goto 79
      y=(y-fmatom(7))/2.556   ! (12)-(7)
      l=5
      goto 66
79    if(ksod.gt.17) goto 65
      y=(y-fmatom(12))/1.362   ! (17)-(12)
      l=6
      goto 66
  65  y=(y-fmatom(17))/1.496   ! (27)-(17)
      l=7
66    do 67 i=2,20
      if(lambda.gt.wmatom(i)) goto 67
      k=i
      goto 68
  67  continue
  68  kk=k-1
      y2=(lambda-wmatom(kk))/(wmatom(k)-wmatom(kk))
      y1=cmatom(kk,l)+y2*(cmatom(k,l)-cmatom(kk,l))
      y2=cmatom(kk,l+1)+y2*(cmatom(k,l+1)-cmatom(kk,l+1))
      mgatom=z1*pg(12)*10.**(y*(y1-y2)-y1-20.)
	dmgatom=mgatom*(dz1/z1+dpg(12))
	ddmgatom=mgatom*ddpg(12)
c.......................................................................

52    	kac=kac+carbon+sodium+mgatom
	kat=kat+dcarbon+dsodium+dmgatom
	kap=kap+ddcarbon+ddsodium+ddmgatom
99999 wavelt=0.			! antes ponia wavelt=lambda
      end
c=======================================================================

      function fkny(i,x)
      real f1(12)/14.47,-169.385,11.65,26.57,31.059,35.31,35.487,
     *5.51,10.36,21.41,37.,25.54/,f2(12)/-2.,21.035,-1.91,-2.9,-3.3,-3.5
     *,-3.6,-1.54,-1.86,-2.6,-3.69,-2.89/
      fkny=f1(i)+x*f2(i)
      if(i-2.lt.0) goto 1
      if(i-2.gt.0) goto 1
      fkny=fkny-0.727*x**2
   1  return
      end
c----------------------------------------------------------------------
	real function um(m,t)
	integer m
	um=1.568399e+5/(t*float(m**2))
	return
	end
c----------------------------------------------------------------------
	real function gff(t,x)
      	gff=1.0828+3.865e-6*t+x*(7.564e-7+(4.92e-10-2.482e-15*t)*t+
     *	x*(5.326e-12+(-3.904e-15+1.879e-20*t)*t))
	return
	end
c-----------------------------------------------------------------------
	real function dum(dum1,t)
c	dum1 es el resultado de calcular um(m,t)
c	dum1= um(m,t)=1.568399e+5/(t*float(m**2))

	dum=-1.*dum1/t
	return
	end
c-----------------------------------------------------------------------

	real function dgff(t,x)
c
c      	gff(t,x)=1.0828+3.865e-6*t+x*(7.564e-7+(4.92e-10-2.482e-15*t)*t+
c     *	x*(5.326e-12+(-3.904e-15+1.879e-20*t)*t))
c
	dgff=3.865e-6+x*(4.92e-10-2.*2.482e-15*t+x*(-3.904e-15+2.*t*
     *  1.879e-20))
c	
	return
	end
c======================================================================


