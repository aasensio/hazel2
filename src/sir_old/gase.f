c	gase calcula las presiones parciales y el peso del ion medio y del neutro medio
c mui=peso medio de los iones
c mun=peso medio de los neutros
c pn=presion de los neutros


      subroutine gase(t,pe,pg,pp,mui,mun,pn)
      parameter (ncontr=28)
      real*4 alfai(ncontr),chi1(ncontr),chi2(ncontr),
     *u0(ncontr),u1(ncontr),u2(ncontr)
      real*4 du0,du1,du2,mui,mun,pesos(ncontr),num_iones,num_neutros,pn
      real*8 pp(*),p(ncontr)
      real*8 t,pe,pg,g1,g2,g3,g4,g5,thetad,g3i
      real*8 a,b,c,d,e,c1,c2,c3,f1,f2,f3,f4,f5,fe,phtot
      real*4 theta,cmol(91),dcmol(91) 


      theta=5040./t
      thetad=5040.d0/t
      call molecb(theta,cmol,dcmol)
      g4=pe*10.**(cmol(1))
      g5=pe*10.**(cmol(2))


c	ahora calculo las funciones de particion u0,u1,u2 y sus derivadas
        do 5 i=1,ncontr
      		iii=i
5     		call neldatb(iii,0.,pesos(i),alfai(i),chi1(i),chi2(i))
6       do 4 i=1,ncontr
      	  iii=i
	  t0=t
     	  call nelfctb(iii,t0,u0(iii),u1(iii),u2(iii),du0,du1,du2)
4       end do
      
        g2=sahadouble(thetad,chi1(1),u0(1),u1(1),pe)   ! p(h+)/p(h)
        g3=1./sahadouble(thetad,0.754,1.,u0(1),pe)     ! p(h-)/p(h)
        pp(10)=g3
        g1=0.
        num_neutros=0.     !numero  neutro/h total (inicializacion para sumar. Falta HI)
        num_iones=0.       !numero iones /h total (inicializacion para sumar. Falta H+ y H-)
        mun=0.             !peso de los neutros
        mui=0.             !peso de los iones
        do i=2,ncontr
          a=sahadouble(thetad,chi1(i),u0(i),u1(i),pe) !a=n(FeII)/n(FeI)
          b=sahadouble(thetad,chi2(i),u1(i),u2(i),pe) !b=n(FeIII)/n(FeII)
	
	  c=1.+a*(1.+b)        !c=n(Fe)/n(FeI)
          p(i)=alfai(i)/c      ! p/ph' for neutral he,li, ... o sea n(FeI)/n(Htot)
	
	  num_neutros=num_neutros+p(i)  !n(FeI)/n(Htot)
	  mun=mun+pesos(i)*p(i)
	
	  num_iones=num_iones+p(i)*a*(1.+b)
	  mui=mui+pesos(i)*p(i)*a*(1.+b) !p(i)*a*(1.+b)=n(FeI)/n(Htot)*n(FeII)/n(FeI)*[1.+n(FeIII)/n(FeII)]
		
	  g1=g1+p(i)*a*(1.+2.*b)   !numero total de e que no vienen de ionizar H
        end do
        pp(2)=p(2)
        pp(3)=p(6)
        pp(4)=p(11)
        pp(5)=p(12)
 
        a=1.+g2+g3
        b=2.*(1.+g2/g5*g4)
        c=g5
        d=g2-g3
        e=g2/g5*g4
        c1=c*b**2.+a*d*b-e*a**2.
        c2=2.*a*e-d*b+a*b*g1	
        c3=-(e+b*g1)
        f1=0.5*c2/c1
        f1=-f1+sign(1.d0,c1)*sqrt(f1**2.-c3/c1)
        f5=(1.-a*f1)/b  !numero de h2/n(h')
        f4=e*f5         !numero de h2+/n(h')
        f3=g3*f1        !numero de h-/n(h')
        f2=g2*f1        !numero de h+/n(h')
        fe=f2-f3+f4+g1  !numero de e/n(h')
	phtot=pe/fe     !kT ne/(ne/nh')=n(h')kt

        if(f5.gt.1.e-4) goto 2
          const6=g5/pe*f1**2
          const7=f2-f3+g1

	  do 3 i=1,5
      	       f5=phtot*const6
	       f4=e*f5
	       fe=const7+f4
	       phtot=pe/fe
3	  continue
      
2       pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe)

        pp(1)=f1   ! p(h)/p(h')
        pp(6)=f2   ! p(h+)/p(h')	
        pp(7)=f5   ! p(h2)/p(h')
        pp(8)=fe   ! pe/p(h')
        pp(9)=pe/(1.38054e-16*t) ! n(e)=pe/kt
		
	num_neutros=num_neutros+f1+f5   !num neutros/n(h')
	mun=(mun+f1+2.*f5)/num_neutros  !peso medio neutros
		
	num_iones=num_iones+f2+f3+f4    !num iones/n(h')
	mui=(mui+f2+f3+2.*f4)/num_iones  !peso medio iones
	
	pn=num_neutros*phtot       !presion de los neutros

      return
      end

