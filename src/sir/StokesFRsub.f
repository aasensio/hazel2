c_______________________________________________________________
c copiada de blends2 (rutina de SIR)
c se basa en la rutina "blends0" asi pues esta disenyado
c igual que la rutina "fstokes0" pero teniendo en cuenta blends
c calcula las funciones respuesta de los perfiles de stokes
c para un angulo heliocentrico de coseno mu. Se supone lte,
c se basa en el programa "zeemanlines" de Wittmann y usa sus rutinas
c modificadas y en el "zeeman" de Rees usando el metodo hermite para integrar
c Basilio Ruiz 23-3-93
c Basilio y Jose Carlos 9-1-95 para pert. aditivas en fi
c Basilio y Jose Carlos 6-2-95 para pert. aditivas en gamma
c Basilio 24-1-2016
c _________________________________________________________________
c
c ist=1 (i); =2 (q); =3 (u); =4 (v)

        subroutine StokesFRsub(stok,rt,rp,rh,rv,rg,rf,rm,rmac,dLPgdPe,dLPgdT)     

        include 'PARAMETER'         !incluye kt,kn,kl,kld

        parameter (kt4=4*kt,kt16=16*kt,kt7=7*kt,kt8=8*kt+2)
        parameter (kl4=4*kl,kld4=4*kld) 

        real*4 melectron,mhidrogeno
        
c new variables to include leelineasii data      
        character*2 atom_all(kl)
        integer istage_all(kl)
        real*4 wlengt_all(kl),zeff_all(kl),energy_all(kl)
        real*4 loggf_all(kl)
        integer mult_all(2,kl)
        character*1 design_all(2,kl)
        real*4 tam_all(2,kl),alfa_all(kl),sigma_all(kl)

c para las funciones respuesta
        real*4 rt(*),rp(*),rh(*),rv(*),rg(*),rf(*),rm(*),dLPgdPe(*),dLPgdT(*)
        real*4 rmac(*)
        real*4 rtlin(kld4),rplin(kld4),rhlin(kld4)
        real*4 rvlin(kld4),rglin(kld4),rflin(kld4)
        real*4 rmlin(kld4)
        real*4 rt4(4,kt),rp4(4,kt),rh4(4,kt)
        real*4 rv4(4,kt),rg4(4,kt),rf4(4,kt)
        real*4 rm4(4,kt)
        real*4 grt(kt),grp(kt),grh(kt),grv(kt)
        real*4 grg(kt),grf(kt),grm(kt)
        real*4 x(kt)

c para los perfiles 
        real*4 atmosmodel(kt8),stok(*),svec(4)

c para la atmosfera
        real*4 t(kt),pe(kt),vtur(kt),h(kt),vz(kt)
        real*4 bp(kt),bt(kt),dbp(kt)
        real*4 gamma(kt),phi(kt),agamma(kt),aphi(kt)
        real*4 vof(kt),tau(kt),taue(kt)
        real*4 voffset,xmu
        real*4 continuoh,conhsra    !,dplnck,dtplanck
        integer*4 mnodos(18),ifiltro

c para la matriz de absorcion y sus derivadas
        real*4 dab(kt16),tk(kt16),pk(kt16)
        real*4 hk(kt16),vk(kt16),gk(kt16)
        real*4 fk(kt16),mk(kt16)
        real*4 dabtot(kt7,kld),gktot(kt7,kld),fktot(kt7,kld)
        real*4 tktot(kt7,kld),pktot(kt7,kld),hktot(kt7,kld)
        real*4 vktot(kt7,kld),mktot(kt7,kld)
        real*4 dabb(16),gkb(16),fkb(16),tkb(16)
        real*4 pkb(16),hkb(16),vkb(16)
        real*4 mkb(16)
        real*4 beta1(kl,kt),beta2(kl,kt)

c para la malla
        real*4 wlengt,wlengt1,lambda,dlamdaii,wvac,wc,c
        integer ist(4),ntau
        integer ntl,nlin(kl),npas(kl),nble(kl)
        real*4 dlamda(kld),dlamda0(kl)
        integer ntls,nlins(kl4),npass(kl4)
        real*4 dlamdas(kld4),dlamda0s(kl4)
c para los parametros atomicos y coeficientes de absorcion
        real*4 loggf,nair,mvdop,meta0,ma
        real*4 y(kt),dyt(kt),dyp(kt),alpha(kt),www
        real*4 table(kt,16)
        real*4 ck5(kt),dk5(kt),ddk5(kt),zeff
        real*4 ckappa,ckappa5,dkappa,dkappa5,ddkappa,ddkappa5
        real*4 piis
        
c para el patron zeeman
        character atom*2
        parameter (mc=20)       !numero maximo de componentes zeeman
        integer mult(2),ji(2),jf(2)
        character design(2)*1           !,multno*6
        real*4 tam(2),abu
        real*4 dlp(mc),dll(mc),dlr(mc),sp(mc),sl(mc),sr(mc)

c para las presiones parciales
        integer ivar(10)
        real*4 pg(99),dpg(99),ddpg(99),pi(10),dpi(10),ddpi(10)
        real*4 pt(kt,10),dpt(kt,10),ddpt(kt,10)
        real*4 pgas(kt),dpgas(kt),ddpgas(kt)      !,ro(kt),ck5_ro(kt)
        real*4 macro

c para la inclusion de RP en RT
        real*4 ax(kt),bx(kt),cx(kt),dx(kt),d1x(kt),d2x(kt),fx(kt)
        real*4 px(kt),qx(kt),rx(kt),sx(kt),tx(kt),wx(kt,kt)
        real*4 kac,kat,kap
        
c para hermite
        real*4 deltae(kt),deltai(kt),delt2i(kt)

        integer :: error_code
	common/Error/error_code
        
c lugares comunes de memoria
        common/Atmosmodel/atmosmodel,ntau                  !se carga en lee_model.f
        common/Malla/ntl,nlin,npas,nble,dlamda             !se carga en lee_malla.f
        common/Malla4/ntls,nlins,npass,dlamdas             !se carga en lee_malla.f
        common/yder/y,dyt,dyp,alpha                        !se carga aqui coef.abs.cont.y su der.t,p 
        common/nlte/nlte                                   ! NLTE
        common/departcoef/beta1,beta2
        common/segunda/tau,taue,deltae,deltai,delt2i       !se carga aqui (para hermite y rnorma)
        common/piis/piis                                   !1./sqrt(3.1415926) se carga aqui (para mvoigt)
        common/offset/voffset                              !se carga en lee_model 
        common/anguloheliocent/xmu                         !se carga en lee_model
        common/numeronodos/mnodos                          !se carga en lee_model
        common/ifiltro/ifiltro
        common/Lineas_all/atom_all,istage_all,wlengt_all,zeff_all,  
     &         energy_all,loggf_all,mult_all,design_all,tam_all,  
     &         alfa_all,sigma_all

        data iprimera/0/
                 
c nble es el numero de componentes de cada linea
        c=2.99792458e+10        !vel. de la luz en cm/seg
        piis=1./sqrt(3.1415926)
        g=xmu*2.7414e+4         !gravedad cm/s^2 en fotosfera solar   
        avog=6.023e23
        epsilon=.95             !cota para corregir RB por presiones Pmag<epsilon*Ptot

        bol=1.3807e-16          !erg/s
        pir=3.1415926
        v0=1e6                  !cm/s
        melectron=9.1094e-24
        mhidrogeno=1.67442e-24
        xmasaproton=1.6526e-24
        avo=6.023e23
        bohr=0.0529177249e-7   !cm
        uma=1.660540e-24
        gas=8.31451e7         !constante de los gases en cgs
        epsilon2=epsilon/(1.-epsilon)  ! Pmag< epsilon2 * Pg

        coc3=1.212121          ! cociente polarizabilidad del H2 con HI
        coc2=.3181818          ! idem para el He


        ntotal=0
        do i=1,ntl
           do j=1,npas(i)
                ntotal=ntotal+1
           end do
        end do
        
        do i=1,4
          ist(i)=1          !sintetizamos los 4 parametros de Stokes
        end do
        
        ists=ist(1)+ist(2)+ist(3)+ist(4)
        ntotal4=ntotal*ists
        
        macro=atmosmodel(8*ntau+1)
        
        if(iprimera.eq.0)then
           do i=1,ntau
              tau(i)=atmosmodel(i)
              taue(i)=10.**(tau(i))
           end do
           do i=1,ntau-1
              deltae(i)=taue(i)-taue(i+1)
           end do 
           do i=2,ntau
              deltai(i)=(tau(i)-tau(i-1))/2.0
              delt2i(i)=deltai(i)*deltai(i)/3.0
           end do

           paso=tau(1)-tau(2)
           iprimera=1
        end if

        x(1)=g*(tau(1)-tau(2))*2.3025851         
        do i=1,ntau-1                             
           x(i+1)=g*(tau(i)-tau(i+1))*2.3025851   
        end do                                   

c leemos la atmosfera
        do i=1,ntau
           t(i)=atmosmodel(i+ntau)
           pe(i)=atmosmodel(i+2*ntau)
           if (pe(i).lt.0.) then
              print*,' '
              print*,'STOP: pe is NEGATIVE at some optical depth.'
              print*,'      Subroutine StokesFRsub: pe(',i,')= ',pe(i)
              print*,' '
              print*,'_______________________________'
              error_code = 1
              return
           end if
           vtur(i)=atmosmodel(i+3*ntau)
           h(i)=atmosmodel(i+4*ntau)
           vz(i) =atmosmodel(i+5*ntau)       !velocidad linea de vision
           vof(i)=vz(i)-voffset
           gamma(i)=atmosmodel(i+6*ntau)
c           print*,'StokesFRsub 182',i,gamma(i)
           phi(i)=atmosmodel(i+7*ntau)
c          agamma(i)=tan(gamma(i)/2.0)
c          aphi(i)=tan(phi(i)/4.0)
           aphi(i)=1. !no qeremos pertb. relativas en fi sino aditivas
           agamma(i)=1. !no qeremos pertb. relativas en gamma sino aditivas
        end do

c calculamos el peso molecular medio pmu
        pmusum=0.0
        asum=0.0
        do i=1,92
           ii=i
           call neldatb(ii,0.,wgt,abu,ei1,ei2)
           pmusum=pmusum+wgt*abu
           asum=asum+abu
        end do

c calculo las presiones parciales (pg) y las derivadas de sus logaritmos
c con respecto a la t (dpg) y, con pe (ddpg)

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

        do i=1,ntau
            ps=pe(i)
            ts=t(i)
            theta=5040./ts
            call gasb(theta,ps,pg,dpg,ddpg)
            if (error_code == 1) return
            do j=1,10
               k=ivar(j)
               pt(i,j)=pg(k)
               dpt(i,j)=dpg(k)
               ddpt(i,j)=ddpg(k)
               pi(j)=pg(k)
               dpi(j)=dpg(k)
               ddpi(j)=ddpg(k)
            end do
            pgas(i)=pg(84)      !presion gaseosa
            dpgas(i)=dpg(84)    !der. log(presion gaseosa) respectp a T
            ddpgas(i)=ddpg(84)  !der. log(presion gaseosa) respectp a pe
            call kappach(5.e-5,ts,ps,pi,dpi,ddpi,ck5(i),dk5(i),ddk5(i))
            cc=avog/pmusum
            kac=ck5(i)*cc
            kat=dk5(i)*cc
            kap=ddk5(i)*cc
            tauk=taue(i)/2./kac/kac

            fx(i)=dpg(84)*pg(84)
            cx(i)=1./(ddpg(84)*pg(84))
            bx(i)=kap/(ddpg(84)*pg(84))
            ax(i)=kat-kap*dpg(84)/ddpg(84)
            d1x(i)=x(i)*tauk   !AQUI BRC
            d2x(i)=x(i)*tauk
            dx(i)=1.+d1x(i)*bx(i)

        end do

        

        rx(ntau)=0. 
        sx(ntau)=0. 
        tx(ntau)=0. 
        do i=1,ntau-1
           rx(i)=(1.- d2x(i+1)*bx(i+1))/dx(i)
           sx(i)=(d2x(i+1)/dx(i))*ax(i+1)
           tx(i)=(d1x(i)/dx(i))*ax(i)
           px(i)=-cx(i)*(tx(i)+fx(i))
        end do

        do i=2,ntau
           qx(i)=-cx(i-1)*(sx(i-1)+tx(i)*rx(i-1))
        end do
        qx(1)=0.

        do i=1,ntau-1
c          deltapei=px(i)*deltat(i)+qx(i+1)*deltat(i+1)
           do j=1,i+1
              wx(i,j)=0.
           end do
           do j=ntau-1,ntau
              wx(i,j)=0.
           end do
           do j=i+2,ntau-2
              r1=1.
              do k=i,j-2
                 r1=r1*rx(k)
              end do
              wxi=-(sx(j-1)+tx(j)*rx(j-1))*r1 
c             deltapei=deltapei+cx(i)*wxi*deltat(j)
              wx(i,j)=cx(i)*wxi
              
           end do
c          deltape(i)=deltapei
        end do


c        call densidad(ntau,tau,t,pe,ck5_ro,pgas,ro) !ck5_ro=kap/ro

c datos de la linea
        ikk0=0
        ikk1=0
        ixx=0
        do 999 iln=1,ntl        !iln=numero de la linea,ntl num.total
           do i=1,npas(iln)
              do k=1,ntau*7
                 dabtot(k,i)=0.
                 tktot(k,i)=0.
                 pktot(k,i)=0.
                 hktot(k,i)=0.
                 gktot(k,i)=0.
                 fktot(k,i)=0.
                 vktot(k,i)=0.
                 mktot(k,i)=0.
              end do
           end do
           do ible=1,nble(iln)   
              ixx=ixx+1
              nxx=nlin(ixx)                    !Including data instead of reading with leelineas
              atom=atom_all(ixx)
              istage=istage_all(ixx)
              wlengt=wlengt_all(ixx)
              zeff=zeff_all(ixx)
              energy=energy_all(ixx)
              loggf=loggf_all(ixx)
              mult(1)=mult_all(1,ixx)
              mult(2)=mult_all(2,ixx)
              design(1)=design_all(1,ixx)
              design(2)=design_all(2,ixx)
              tam(1)=tam_all(1,ixx)
              tam(2)=tam_all(2,ixx)
              alfa=alfa_all(ixx)
              sigma=sigma_all(ixx)
         
c              if(nxx.eq.0)then
c                 nxx=nlin(ixx-1)
c                 call leelineasii(nxx,atom,istage,wlengt,zeff,energy,
c     &                         loggf,mult,design,tam,alfa,sigma)
c                 loggf=-20.
c                 wlengt=5000.
c              else
c                 call leelineasii(nxx,atom,istage,wlengt,zeff,energy,
c     &                         loggf,mult,design,tam,alfa,sigma)
c              end if   
              gf=1.e1**(loggf)
              if(ible.eq.1)wlengt1=wlengt
              dlamda0(iln)=wlengt
              continuoh=conhsra(wlengt1)

c parametros atomicos
c llamo a atmdat que devuelve weight (peso molecular),abu (abunda
c cia), chi1,chi2 (pot.de ionizacion del atomo neutro e ion),
c u1,u2,u3 (funciones de particion atomo neutro,ion,ion2).
        call atmdatb(atom,0.,nel,weight,abu,chi10,chi20,u1,u2,u3,
     &               du1,du2,du3)

c calculo la l.d.o. en el vacio (wvac), con la function refrax
        nair=1.0004
        if(wlengt.ge.1800.)nair=refrax(wlengt*1.e-4,15.,760.,0.)
        lambda=wlengt*1.e-8     !l.d.o. en el aire en cm.
        wvac=nair*lambda        !l.d.o. en el vacio en cm.
        wc=wvac/c               !inverso de la frecuencia (segundos)

c calculo de terminos constantes.
        weinv=1./weight
        croot=1.66286e+8*weinv  !2r/m  para anchura doppler(cm**2/s**2)
        dlo=4.6686e-5*wvac*wvac !separac.de l.d.o=dlo*h*(m1*g1-m2*g2)
        crad=.22233/wvac        !(l.d.o.*coef.clasico ensanch. natural)
        eta00=1.49736e-2*gf*abu*wvac !para el coef. de absor.de la linea

c calculo chydro (gamma6)coef.van der waals ensanch.colisional
c corregido con zeff (termino 'semiempirico' caca de la vaca)
c (se supone despreciable el debido a stark)
c si el damping asi calculado es mayor que 3 lo reescalaremos

        if(alfa.eq.0.or.sigma.eq.0)then

           chi1=chi10
           if(istage.eq.2)chi1=chi20
           eupper=energy+1.2398539e-4/wvac      
           ediff1=amax1(chi1-eupper-chi10*float(istage-1),1.0)

           ediff2=amax1(chi1-energy-chi10*float(istage-1),3.0)
           chydro=lambda*1.e1**(+.4*alog10(1./ediff1**2-1./ediff2**2)-
     &            12.213)*5.34784e+3

           chydro=chydro*zeff

           if(istage.eq.2)chydro=chydro*1.741

        else

c xmu1 es la masa reducida del hidrogeno y el atomo.Los indices 2 y 3 hacen
c referencia al helio neutro y al hidrogeno molecular.
   
           xmu1=uma*(1.008*weight)/(1.008+weight)
           xmu2=uma*(4.0026*weight)/(4.0026+weight)
           xmu3=uma*(2.016*weight)/(2.016+weight)

           arr=2.-alfa*.5-1.
           gammaf=1.+(-.5748646+(.9512363+(-.6998588+(.4245549-
     &              .1010678*arr)*arr)*arr)*arr)*arr
           vv=(1.-alfa)/2.

           beta=lambda*2*(4./pir)**(alfa/2.)*gammaf*(v0**alfa)*sigma*
     &          ((8.*bol/pir)**vv)

        endif

c subniveles zeeman
c calculo la parte entera "ji" y la parte fraccionaria "jf" del
c momento angular total (tam), necesarios para zeeman
        do i=1,2
           ji(i)=int(tam(i))
           jf(i)=int(10*(tam(i)-ji(i)))
        end do

c llamo a zeeman para que calcule el numero de transiciones pi,
c l,r (np,nl,nr),desplazamientos de cada componente (dlp,dll,dlr)
c e intensidad (sp,sl,sr).dlo se calculo en el apartado @5

        if(wlengt.gt.15652.7.and.wlengt.lt.15652.9)then
c            !es la linea infrarroja con nivel superior en acoplamiento jk
             !hay que indicar los factores de lande de los dos niveles.
             !el J de cada nivel ya lo coje del fichero LINEAS
c            print*,'entro en zeeman_jk'

c            open(66,file='gfactors')
c            read(66,*) xgg1,xgg2
c            close(66)

c            call zeeman_jk(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
c     &      sp,sl,sr,1.510,1.499)

             call zeeman_jk(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
     &       sp,sl,sr,1.45,1.45)
        else
             call zeeman(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
     &       sp,sl,sr)
        endif

c estratificacion en tau 5000
c genero una tabla lineal en logaritmo de tau a l.d.o.=5000.angs.
c desde tauini a taufin, con ntau valores.

        do 71 i=1,ntau
            ps=pe(i)
            ts=t(i)
            theta=5040./ts

            do j=1,10
               k=ivar(j)
               pi(j)=pt(i,j)
               dpi(j)=dpt(i,j)
               ddpi(j)=ddpt(i,j)
               pg(k)=pi(j)
               dpg(k)=dpi(j)
               ddpg(k)=ddpi(j)
            end do

c calculo el coeficiente de absorcion del continuo por cm**3
c (ckappa) y sus derivadas con respecto a t y pe:dkappa,
c ddkappa y lo divido por dicho coef. evaluado a 5000.a
c calculo el ckappa para lambda=(5000 a=5.e-5 cm),(ckappa5)

            call kappach(lambda,ts,ps,pi,dpi,ddpi,ckappa,
     &                       dkappa,ddkappa)
            ckappa5=ck5(i)
            dkappa5=dk5(i)
            ddkappa5=ddk5(i)

            ckappa=ckappa/ckappa5
            dkappa=(dkappa-ckappa*dkappa5)/ckappa5
            ddkappa=(ddkappa-ckappa*ddkappa5)/ckappa5
            dkappa5=dkappa5/ckappa5
            ddkappa5=ddkappa5/ckappa5

c       calculo vdop=(delta doppler*c/l.d.o.) y su derivada
            vdop=sqrt(croot*ts+vtur(i)**2) !cm/seg.
            vdop2=sqrt((2.*gas*ts)/weight+vtur(i)**2)
            dvdop=croot/vdop/2.
            dvdop2=gas/(vdop2*weight)
            mvdop=vtur(i)/vdop

c calculo el coeficiente absorcion linea en cada tau y sus derivadas
c calculo eta0=(coef.absorcion linea/coef. absorcion continuo)
c necesito calcular la fraccion de atomos del elemento en el
c el nivel de la transicion respecto al numero total de atomos
c del elemento en cualquier estado de ionizacion.asi, tengo que
c utilizar las ecuacione de saha y boltzmann.
c llamare u12 al cociente de poblaciones entre el estado de ioniza
c cion 1 (neutro) y 2.igualmente u23 entre los iones 2 y 3
c necesito llamar a nelfctb (definida en atmdatb con un 'entry')
c para calcular las funciones de particion y sus derivadas a cada
c temperatura
           call nelfctb(nel,ts,u1,u2,u3,du1,du2,du3)

c corrijo los potenciales de ionizacion chi10 y chi20 con
c un termino proporcional a la raiz cubica de la densidad de e-

           rcu=(pg(91))**(1./3.)
           chi1=chi10-6.96e-7*rcu
           chi2=chi20-1.1048e-6*rcu

           u12=saha(theta,chi1,u1,u2,ps)      !n2/n1
           du12=u12*dsaha(theta,chi1,du1,du2) !derivada de u12 con t
           ddu12=-1.*u12/ps        !    "    "   "  con pe
           u23=saha(theta,chi2,u2,u3,ps)      !n3/n2
           du23=u23*dsaha(theta,chi2,du2,du3) !derivada de u23 con t
           ddu23=-1.*u23/ps                !    "    "   "  con pe
           u33=1.+u12*(1.+u23)                !(n1+n2+n3)/n1
           du33=du12*(1.+u23)+u12*du23
           ddu33=ddu12*(1.+u23)+u12*ddu23

           eta0=eta00*1.e1**(-theta*energy)/(u1*vdop*u33*ckappa5)
           deta0=alog(10.)*theta/ts*energy-du1-dvdop/vdop-du33/u33
           deta0=eta0*(deta0-dkappa5)
           ddeta0=eta0*(-ddu33/u33)
           ddeta0=ddeta0-eta0*ddkappa5

           if(istage.eq.2)then
                eta0=eta0*u1/u2*u12
                deta0=deta0*u1/u2*u12+eta0*(du1-du2+du12/u12)
                ddeta0=ddeta0*u1/u2*u12+eta0*(ddu12/u12)
           end if

c introduzco la correccion por emision estimulada
           corre=1.-exp(-1.4388/(ts*wvac))
           dcorre=(corre-1.)*(1.4388/(ts*ts*wvac))
           deta0=deta0*corre+eta0*dcorre
           ddeta0=ddeta0*corre
           eta0=eta0*corre
           meta0=-eta0*mvdop/vdop

c calculo el damping 'a' mediante WITTMANN.
c pg(1)=p(h)/p(h'),pg(90)=pe/p(h'),pg(91)=densidad de e-
c pg(2)=p(he)/p(h'),pg(89)=p(h2)/p(h')
c   a=(chydro*(pg(1)/pg(90)*pg(91))*t(i)**0.3*((.992093+weinv)**.3
c  &  +.6325*pg(2)/pg(1)*(.2498376+weinv)**.3+.48485*pg(89)/pg(1)*
c  &  (.4960465+weinv)**.3)+crad)/(12.5663706*vdop)


        if(sigma.eq.0.or.alfa.eq.0)then

           if(pg(1).gt.1.e-20)then !!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0
        
              aj=chydro*(pg(1)/pg(90)*pg(91))*ts**0.3
              ai=(.992093+weinv)**.3+.6325*pg(2)/pg(1)*(.2498376+weinv)
     &           **.3+.48485*pg(89)/pg(1)*(.4960465+weinv)**.3
              a=(aj*ai+crad)/(12.5663706*vdop)
              ma=-a*mvdop/vdop
        
              daj=aj*(dpg(1)-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
              ddaj=aj*(ddpg(1)-ddpg(90)+ddpg(91)) !derivada de aj con p

              dai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(dpg(2)-dpg(1))
     &            +.48485*pg(89)/pg(1)*(.4960465+weinv)**.3*(dpg(89)
     &             -dpg(1))
              ddai=.6325*pg(2)/pg(1)*(.2498376+weinv)**.3*(ddpg(2)
     &             -ddpg(1))+.48485*pg(89)/pg(1)*(.4960465+weinv)**
     &             .3*(ddpg(89)-ddpg(1))

              da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
              dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)

           else !!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0

              aj=chydro*(1./pg(90)*pg(91))*ts**0.3
              ai=.6325*pg(2)*(.2498376+weinv)**.3+
     &      .48485*pg(89)*(.4960465+weinv)**.3
              a=(aj*ai+crad)/(12.5663706*vdop)
              ma=-a*mvdop/vdop
        
              daj=aj*(-dpg(90)+dpg(91)+.3/ts) !derivada de aj cont
              ddaj=aj*(-ddpg(90)+ddpg(91)) !derivada de aj con p
              dai=.6325*pg(2)*(.2498376+weinv)**.3*dpg(2)+
     &      .48485*pg(89)*(.4960465+weinv)**.3*dpg(89)
              ddai=.6325*pg(2)*(.2498376+weinv)**.3*ddpg(2)+
     &       .48485*pg(89)*(.4960465+weinv)**.3*ddpg(89)

c              AQUI!!! CREO QUE FALTAN ESTAS DOS SENTENCIAS
	    da=(ai*daj+aj*dai)/(12.5663706*vdop)-a*dvdop/vdop
	    dda=(ai*ddaj+aj*ddai)/(12.5663706*vdop)
c              HASTA AQUI!!! 

           end if!!!!!!!! CORRECCION PARA EVITAR DAMPING NULO POR PG(1)=0

       else

c       calculo del damping 'a' mediante BARKLEM.       
c       d son las derivadas totales respecto a la temperatura.
c       dd son las derivadas totales respecto a la presion.

           dam=beta*(pg(91)/pg(90))*(pg(1)*xmu1**(-vv)+coc2*pg(2)*
     &        (xmu2**(-vv))+coc3*pg(89)*xmu3**(-vv))*ts**vv
        
           a=(1./(4.*pir))*(crad/vdop+dam/vdop2) 
           ma=-a*mvdop/vdop         

           ddam=dam*(dpg(91)-dpg(90))+beta*(pg(91)/pg(90))*(xmu1**(-vv)*
     &       pg(1)*dpg(1)+coc2*xmu2**(-vv)*pg(2)*dpg(2)+coc3*pg(89)*
     &       dpg(89)*xmu3**(-vv))*ts**(vv)+vv*(dam/ts)

           dddam=dam*(ddpg(91)-ddpg(90))+beta*(pg(91)/pg(90))*(xmu1**
     &           (-vv)*pg(1)*ddpg(1)+coc2*xmu2**(-vv)*pg(2)*ddpg(2)+
     &           coc3*pg(89)*ddpg(89)*xmu3**(-vv))

           da=(crad/(4.*pir))*(-dvdop/(vdop**2.))+(1./(4.*pir))*(ddam*
     &        vdop2-dam*dvdop2)/(vdop2**2.)

           dda=(1./(4.*pir))*(dddam/vdop2)

       endif

c calculo la funcion de planck en lambda y su derivada con t
c la function "dtplanck" calcula la derivada de la f. de planck
c con la temperatura.
                www=wlengt1*1.e-8


c Add the effect of departure coefficients. If the line is in LTE (the usual case)
c they will be equal to 1
                blow=beta1(ixx,i)
                bratio=blow/beta2(ixx,i)
                
                call planck2(t(i),www,bratio,bp(i),bt(i))

                eta0=eta0*blow
                deta0=deta0*blow
                ddeta0=ddeta0*blow
                meta0=meta0*blow

c                bp(i)=dplnck(t(i),www)  !cuerpo negro
c                bt(i)=dtplanck(t(i),www)               
                y(i)=ckappa              !ckappa/ckappa5
                dyt(i)=dkappa
                dyp(i)=ddkappa
                table(i,1)=dkappa        !derivada con t
                table(i,2)=ddkappa       !    "     "  pe
                table(i,3)=eta0    !kap. linea/kcont.5000
                table(i,4)=deta0   !derivada con t
                table(i,5)=ddeta0  !derivada con pe
                table(i,6)=wc*vdop !despl. doppler en l.d.o(cm)
                table(i,7)=vz(i)/vdop !velocidad eje z (u.doppler)
                table(i,8)=dkappa5       !derivada con t de kappacon5
                table(i,9)=ddkappa5      !  "       "  pe      "
                table(i,10)=a      !damping
                table(i,11)=da   !derivada del damping con t
                table(i,12)=dda    !derivada del damping con p
                table(i,13)=dvdop/vdop   !derivada de log(vdop) con t
                table(i,14)=mvdop/vdop   !   "      "    "      con mic
                table(i,15)=meta0        !derivada de eta0 con la micro
                table(i,16)=ma           !derivada de a    con la micro
71      continue            !do en log(tau)

        call derivacuad(bp,dbp,ntau,tau)

        amaxim=0.
        do i=1,ntau
           if(table(i,10).ge.amaxim)amaxim=table(i,10)
        end do

        if(amaxim.gt.3.)then
           do i=1,ntau
              table(i,10)=(3.*table(i,10))/amaxim
              table(i,11)=(3.*table(i,11))/amaxim
              table(i,12)=(3.*table(i,12))/amaxim
              table(i,16)=(3.*table(i,16))/amaxim
           end do
        end if

c muestreo en lambda , calculo las l.d.o. de cada punto
        do 10 i=1,npas(iln)
           ikk=ikk0+i
           if(ible.eq.nble(iln).and.i.eq.npas(iln))ikk0=ikk0+npas(iln)
           dlamdaii=(dlamda(ikk)+(wlengt1-wlengt)*1.e3)*1.e-11  !en cm.

c calculo etar,etal,etap y sus derivadas
           do 101 j=1,ntau     !do en tau
              k=(j-1)*7+1
              a=table(j,10)     !damping
              dldop=table(j,6) !desplazamiento doppler en cm
              dvcam=table(j,7) !campo de velocidades unidades doppler
              v=dlamdaii/dldop-dvcam !(l.d.o.+c.vel.) en unidades doppler
              hh=h(j)/dldop
              if(ible.eq.1)then
                t0=y(j)         !coef. absor. continuo/c.a.c 5000
                t1=table(j,1)     !derivada de t0 con t
                t2=table(j,2)     !  "         "      pe
              else
                t0=0.
                t1=0.
                t2=0.
              end if
              t3=table(j,3)     !coef. absor. linea/c.a.c 5000
              t4=table(j,4)     !derivada de t3 con t
              t5=table(j,5)     !  "         "   "  pe
              t15=table(j,15)   !  "         "   "  micro
              t11=table(j,11)   !derivada del damping con t
              t12=table(j,12)   !  "       "     "     "  pe
              t16=table(j,16)   !  "       "     "     "  micro
              t13=table(j,13)   !derivada de log(vdop) con t
              t14=table(j,14)   !  "      "      "      "  micro
              sg=sin(gamma(j))
              cg=cos(gamma(j))
              s2g=sg*sg
              c2g=1.+cg*cg
              scg=sg*cg
              sf=sin(2.*phi(j))
              cf=cos(2.*phi(j))
        
c              print*,'tau(',j,')=',taue(j),'  damp =  ',a,'  v = ',v
              call mvoigt(nr,dlr,sr,a,v,hh,t13,t14,wc,dldop,etar,vetar,
     &              getar,ettar,ettvr,esar,vesar,gesar,essar,essvr,
     &              ettmr,essmr)
              call mvoigt(nl,dll,sl,a,v,hh,t13,t14,wc,dldop,etal,vetal,
     &              getal,ettal,ettvl,esal,vesal,gesal,essal,essvl,
     &              ettml,essml)
              call mvoigt(np,dlp,sp,a,v,hh,t13,t14,wc,dldop,etap,vetap,
     &              getap,ettap,ettvp,esap,vesap,gesap,essap,essvp,
     &              ettmp,essmp)

c calculamos la matriz de absorcion
        
              tm=.5*(etar+etal)
              tn=.5*(etar-etal)
              sm=.5*(esar+esal)
              sn=.5*(esar-esal)
              tpm=.5*(etap-tm)
              spm=.5*(esap-sm)

              fi=t0 + t3 * .5 * ( etap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,dabb)
              do iii=1,7
                 jjj=k+iii-1
                 dabtot(jjj,i)=dabtot(jjj,i)+dabb(iii)
              end do
c calculamos la derivada de la matriz de absorcion con gamma
              if(mnodos(6).ne.0)then

              fi=  2.* t3 * tpm * scg
              fq= fi * cf
              fu= fi * sf
              fv=-t3 * tn * sg

              fi1= 2.* t3 * spm * scg
              fq1=fi1 * cf 
              fu1=fi1 * sf 
              fv1=-t3 * sn * sg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,gkb)
              do iii=1,7
                 jjj=k+iii-1
                 gktot(jjj,i)=gktot(jjj,i)+gkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con phi
              if(mnodos(7).ne.0)then

              fi=  0.
              fq= -2.* t3 * tpm * s2g * sf
              fu=  2.* t3 * tpm * s2g * cf
              fv=  0.

              fq1=-2.* t3 * spm * s2g * sf 
              fu1= 2.* t3 * spm * s2g * cf 
              fv1= 0.

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,fkb)
              do iii=1,7
                 jjj=k+iii-1
                 fktot(jjj,i)=fktot(jjj,i)+fkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con t
              if(mnodos(1).ne.0)then

              detar=t4*etar+t3*(ettar*t11+ettvr) !cambio el signo +ettvr
              detal=t4*etal+t3*(ettal*t11+ettvl)
              detap=t4*etap+t3*(ettap*t11+ettvp)
              desar=t4*esar+t3*(essar*t11+essvr)
              desal=t4*esal+t3*(essal*t11+essvl)
              desap=t4*esap+t3*(essap*t11+essvp)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=t1 + .5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,tkb)
              do iii=1,7
                 jjj=k+iii-1
                 tktot(jjj,i)=tktot(jjj,i)+tkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con mic
              if(mnodos(3).ne.0)then

              detar=t15*etar+t3*(ettar*t16+ettmr)!cambio el signo +ettmr
              detal=t15*etal+t3*(ettal*t16+ettml)
              detap=t15*etap+t3*(ettap*t16+ettmp)
              desar=t15*esar+t3*(essar*t16+essmr)
              desal=t15*esal+t3*(essal*t16+essml)
              desap=t15*esap+t3*(essap*t16+essmp)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=.5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,mkb)
              do iii=1,7
                 jjj=k+iii-1
                 mktot(jjj,i)=mktot(jjj,i)+mkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con pe
              if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &        mnodos(6).ne.0)then

              detar=t5*etar+t3*(ettar*t12) 
              detal=t5*etal+t3*(ettal*t12)
              detap=t5*etap+t3*(ettap*t12)
              desar=t5*esar+t3*(essar*t12)
              desal=t5*esal+t3*(essal*t12)
              desap=t5*esap+t3*(essap*t12)
              tm=.5*(detar+detal)
              tn=.5*(detar-detal)
              sm=.5*(desar+desal)
              sn=.5*(desar-desal)
              tpm=.5*(detap-tm)
              spm=.5*(desap-sm)

              fi=t2 + .5 * ( detap*s2g + tm*c2g )
              fq=tpm * s2g * cf
              fu=tpm * s2g * sf
              fv=tn * cg

              fq1=spm * s2g * cf 
              fu1=spm * s2g * sf 
              fv1=sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,pkb)
              do iii=1,7
                 jjj=k+iii-1
                 pktot(jjj,i)=pktot(jjj,i)+pkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con el campo
              if(mnodos(4).ne.0)then

              tm=.5*(getar+getal)
              tn=.5*(getar-getal)
              sm=.5*(gesar+gesal)
              sn=.5*(gesar-gesal)
              tpm=.5*(getap-tm)
              spm=.5*(gesap-sm)

              fi=t3 * .5 * ( getap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,hkb)
              do iii=1,7
                 jjj=k+iii-1
                 hktot(jjj,i)=hktot(jjj,i)+hkb(iii)
              end do
              end if

c calculamos la derivada de la matriz de absorcion con la velocidad
              if(mnodos(5).ne.0)then

              tm=.5*(vetar+vetal)
              tn=.5*(vetar-vetal)
              sm=.5*(vesar+vesal)
              sn=.5*(vesar-vesal)
              tpm=.5*(vetap-tm)
              spm=.5*(vesap-sm)

              fi=t3 * .5 * ( vetap*s2g + tm*c2g )
              fq=t3 * tpm * s2g * cf
              fu=t3 * tpm * s2g * sf
              fv=t3 * tn * cg

              fq1=t3 * spm * s2g * cf 
              fu1=t3 * spm * s2g * sf 
              fv1=t3 * sn * cg 

              call matabs(fi,fq,fu,fv,fq1,fu1,fv1,vkb)
              do iii=1,7
                 jjj=k+iii-1
                 vktot(jjj,i)=vktot(jjj,i)+vkb(iii)
              end do
              end if
        
101          continue !fin do en tau(estamos aun dentro del do en lamda)
10        continue      !fin del do en lambda(pasos)
        end do   !fin del do en blends


        do 9 i=1,npas(iln)
           call matabs2(dabtot,i,ntau,dab)
           if(mnodos(1).ne.0)call matabs2(tktot,i,ntau,tk)
           if(mnodos(1).ne.0.or.mnodos(2).ne.0.or.mnodos(4).ne.0.or.
     &      mnodos(6).ne.0)call matabs2(pktot,i,ntau,pk)
           if(mnodos(3).ne.0)call matabs2(mktot,i,ntau,mk)
           if(mnodos(4).ne.0)call matabs2(hktot,i,ntau,hk)
           if(mnodos(5).ne.0)call matabs2(vktot,i,ntau,vk)
           if(mnodos(6).ne.0)call matabs2(gktot,i,ntau,gk)
           if(mnodos(7).ne.0)call matabs2(fktot,i,ntau,fk)

c       call delospl2(ntau,tau,dab,tk,pk,hk,vk,gk,fk,mk,bp,bt,svec,
c     &  rt4,rp4,rh4,rv4,rg4,rf4,rm4,npaso,ntaun,taun,tauen,bpn,dbdt,s0)

c          call lin(bp,dbp,dab,ntau,svec,kt,bt,tk,pk,hk,
c     &              vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos)
           call hermite(bp,dbp,dab,ntau,svec,kt,bt,tk,pk,hk,
     &              vk,gk,fk,mk,rt4,rp4,rh4,rv4,rg4,rf4,rm4,mnodos)

        if (error_code == 1) return

c introducimos el paso en tau para pasar la integral sobre las f. resp.
c a sumatorio y normalizmos por el continuo


        ikk1=ikk1+1
        ikk4=ikk1-ntotal
        ikk5=ikk1-ntotal
        do isv=1,4
           ikk5=ikk5+ntotal
           if(ist(isv).eq.1)then
              ikk4=ikk4+ntotal

              if(mnodos(2).ne.0.or.mnodos(1).ne.0.or.mnodos(4).ne.0.or.
     &          mnodos(6).ne.0)then
                 do kk=1,ntau
                    grp(kk)=rp4(isv,kk)
                 end do 
              end if

              if(mnodos(1).ne.0)then
                 do kk=1,ntau   !introduccion de grp en grt
                    if(mnodos(2).eq.0)then
                       suma=0.
                       do kj=1,kk-2
                          suma=suma+wx(kj,kk)*grp(kj)
                       end do
                       correc=grp(kk)*px(kk)+suma
                       if(kk.gt.1)correc=correc+grp(kk-1)*qx(kk)
                       grt(kk)=rt4(isv,kk)+correc
                    else
                       grt(kk)=rt4(isv,kk)
                    end if
                 end do 
                 call rnorma(ntau,continuoh,grt)         
c                 call nodos(grt,ntau,tau,t,mnodos(1)) 
              end if
     
              if(mnodos(2).ne.0)then
                 call rnorma(ntau,continuoh,grp)
c                 call nodosp(grp,ntau,tau,pe,mnodos(2))  !OJO nodosp escala conla presion en el ultimo nodo
              end if

              if(mnodos(3).ne.0)then
                 do kk=1,ntau
                    grm(kk)=rm4(isv,kk)
                 end do 
                 call rnorma(ntau,continuoh,grm)
c                 call nodos(grm,ntau,tau,vtur,mnodos(3))   
              end if
              if(mnodos(4).ne.0)then

c                 if(mnodos(2).eq.0.and.ipgmag.eq.1)then
c                   do kk=1,ntau
c                      correc=sin(gamma(kk))
c                      correc=correc*correc*h(kk)/4./3.1415926  !B*sin^2(g)/4./pi
c                      corrijorb=0.
c                      if(correc*h(kk)/2..lt.epsilon2*pgas(kk))corrijorb=1.
c                      correc=correc*corrijorb*cx(kk)
c                      grh(kk)=rh4(isv,kk)-correc*rp4(isv,kk)
c                  end do 
c                 else
                   do kk=1,ntau
                      grh(kk)=rh4(isv,kk)
                   end do 
c                 end if

                 call rnorma(ntau,continuoh,grh)
c                 call nodos(grh,ntau,tau,h,mnodos(4))
              end if
              if(mnodos(5).ne.0)then
                 do kk=1,ntau
                    grv(kk)=rv4(isv,kk)
                 end do 
                 call rnorma(ntau,continuoh,grv)
c                 call nodos(grv,ntau,tau,vof,mnodos(5))
              end if
              if(mnodos(6).ne.0)then
                   do kk=1,ntau
c                     grg(kk)=rg4(isv,kk)*2./(1.0+agamma(kk)*agamma(kk))
                      grg(kk)=rg4(isv,kk)
                   end do 
c                 end if
                 call rnorma(ntau,continuoh,grg)
c                 call nodos(grg,ntau,tau,agamma,mnodos(6))
              end if
              if(mnodos(7).ne.0)then
                 do kk=1,ntau
c                   grf(kk)=rf4(isv,kk)*4./(1.+aphi(kk)*aphi(kk))
                    grf(kk)=rf4(isv,kk)
                 end do
                 call rnorma(ntau,continuoh,grf)
c                 call nodos(grf,ntau,tau,aphi,mnodos(7))
              end if


c las f. respuesta salen ordenadas en longitud de onda,perfil y tau
c rt(tau1:i(l1,l2,...),q(l1,..),....v(l1....);tau2:i......)
              do j=1,mnodos(1)
                 iktt=(j-1)*ntotal4+ikk4
                 rt(iktt)=grt(j)
              end do   
              do j=1,mnodos(2)
                 iktt=(j-1)*ntotal4+ikk4
c Como I=I(T,Pe) -> dI/dT = (dI/dT)_Pe * dT + (dI/dPe)_T * dPe
c El signo menos no tengo claro de donde sale
                 if (j == mnodos(2)) then
                     rp(iktt) = 0.0
                  else
                     rp(iktt)=grp(j) 
                  endif
                  dLPgdT(j) = dpgas(j)
                  dLPgdPe(j) = ddpgas(j)
              end do

              do j=1,mnodos(3)
                 iktt=(j-1)*ntotal4+ikk4
                 rm(iktt)=grm(j)
              end do   
              do j=1,mnodos(4)
                 iktt=(j-1)*ntotal4+ikk4 !ojo las perturbaciones son relativas a
                 rh(iktt)=grh(j)         !los parametros en z no a linea vision
              end do   
              do j=1,mnodos(5)
                 iktt=(j-1)*ntotal4+ikk4
                 rv(iktt)=grv(j)
              end do   
              do j=1,mnodos(6)
                 iktt=(j-1)*ntotal4+ikk4
                 rg(iktt)=grg(j)
              end do   
              do j=1,mnodos(7)
                 iktt=(j-1)*ntotal4+ikk4
                 rf(iktt)=grf(j)
              end do  
 
              stok(ikk4)=svec(isv)/continuoh
              rmac(ikk4)=svec(isv)/continuoh !lo copiamos para convolucionarlo con la derivada de la macro
              end if
            end do
9         continue      !fin del do en lambda
999     continue        !fin del do en lineas


c convolucionamos con la macro (y la PSF si existe)
        if (ifiltro .eq. 1 .or. macro .gt. 0)then
           k1=0
	   do i=1,4
      	      do klin=1,ntl
	         k1=k1+1
	         dlamda0s(k1)=dlamda0(klin)
	      end do
	   end do
	      
         !   call deconv(stok,1,ntls,npass,dlamda0s,dlamdas,macro)
         !   call deconv2(rmac,1,ntls,npass,dlamda0s,dlamdas,macro)

         !   call deconRF(mnodos(1),rt,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(2),rp,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(3),rm,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(4),rh,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(5),rv,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(6),rg,dlamda0s,ntotal4,macro)
         !   call deconRF(mnodos(7),rf,dlamda0s,ntotal4,macro)
        end if
        
        return
        end
c __________________________________________________________________________
        subroutine deconRF(mnod,RF,dlamda0s,ntotal4,macro) 
        include 'PARAMETER'         !incluye kt,kn,kl,kld

        parameter (kl4=4*kl,kld4=4*kld) 
        integer ntls,nlins(kl4),npass(kl4)
        real*4 dlamda0s(*),dlamdas(kld4),macro,RF(*)
        real*4 RFlin(kld4)
        
        common/Malla4/ntls,nlins,npass,dlamdas             !se carga en lee_malla.f
        
        do i=1,mnod
          k1=(i-1)*ntotal4
          k2=k1+ntotal4
          do j=k1,k2
            RFlin(j-k1+1)=RF(j)
          end do  
          call deconv(RFlin,1,ntls,npass,dlamda0s,dlamdas,macro)
          do j=k1,k2
            RF(j)=RFlin(j-k1+1)
          end do
        end do   
        
        return
        end

c __________________________________________________________________________
c matabs rutina que llena la matriz de absorcion
        subroutine matabs(fi,fq,fu,fv,fq1,fu1,fv1,dab)

        real*4 dab(*)

        dab(1)=fi
        dab(2)=fq
        dab(3)=fu
        dab(4)=fv

c       dab(5)=fq
c       dab(6)=fi
c       dab(7)=-fv1
c       dab(8)=fu1

c       dab(9)=fu
c       dab(10)=fv1
c       dab(11)=fi
c       dab(12)=-fq1

c       dab(13)=fv
c       dab(14)=-fu1
c       dab(15)=fq1
c       dab(16)=fi

        dab(5)=fu1
        dab(6)=fq1
        dab(7)=fv1

        return
        end
c__________________________________________________________________________
c matabs2 rutina que construye la matriz de absorcion a partir de los 7
        subroutine matabs2(dab7,ifrec,ntau,dab)

        include 'PARAMETER'  !por kt y kld
        parameter (kt7=7*kt)

        real*4 dab7(kt7,kld)
        real*4 dab(*)

        do jj=1,ntau
           j=16*(jj-1)
           j7=7*(jj-1)

           fi=dab7(j7+1,ifrec)
           fq=dab7(j7+2,ifrec)
           fu=dab7(j7+3,ifrec)
           fv=dab7(j7+4,ifrec)
           fu1=dab7(j7+5,ifrec)
           fq1=dab7(j7+6,ifrec)
           fv1=dab7(j7+7,ifrec)

           dab(j+1)=fi
           dab(j+2)=fq
           dab(j+3)=fu 
           dab(j+4)=fv

           dab(j+5)=fq
           dab(j+6)=fi
           dab(j+7)=-fv1
           dab(j+8)=fu1 

           dab(j+9)=fu
           dab(j+10)=fv1
           dab(j+11)=fi
           dab(j+12)=-fq1

           dab(j+13)=fv
           dab(j+14)=-fu1
           dab(j+15)=fq1
           dab(j+16)=fi
        end do

        return
        end
c _________________________________________________________________