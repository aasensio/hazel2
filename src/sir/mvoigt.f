	subroutine MVOIGT(nr,dlr,sr,a,v,hh,t13,t14,wc,dldop,etar,vetar,
     &      getar,ettar,ettvr,esar,vesar,gesar,essar,essvr,ettmr,essmr)

	REAL*4 WC
	REAL DLR(20),SR(20)
	common/piis/piis


c dlr es el desplazamiento en cm/gauss de la componete r Zeemann
c dlr es negativo luego debo cambiar de signo la exprsion
	      etar=0.    !eta para luz circular derecha
	      vetar=0.   !parcial de eta con el campo de velocidad(s/cm)
	      getar=0.   !parcial de eta con el campo magnetico.
	      ETTAR=0.   !parcial de eta con a
	      ETTVR=0.   !parcial de eta con v
	      ettmr=0.   !parcial de eta con mic
	      ESAR=0.    !perfil anomalo
	      VESAR=0.   !su derivada con la velocidad
	      GESAR=0.   !su derivada con el campo
	      ESSVR=0.   !parcial de esar con v
	      ESSAR=0.   !parcial de esar con a
	      essmr=0.   !parcial de esar con mic

	      w1=wc/dldop
		
  	      do 102 ir=1,nr   !do en el numero de componentes r zeeman
  		 ver=v+dlr(ir)*HH
c	         H=VOIGT(a,ver,0)     !funcion de Voigt
c	         F=VOIGT(a,ver,1)     !funcion de Faraday
	         call voigt2(a,ver,H,F)

	         HV=-2.*ver*H+4.*a*F  !derivada con v de H
	         FV=piis-a*H-2.*ver*F  !derivada con v de F
                 FA=hv/2.             !derivada con a de F
		 HA=-2.*fv    !derivada con a de H
		 etar=etar+sr(ir)*H
	         esar=esar+sr(ir)*F
	         ettaR=ettaR+sr(ir)*HA
	         essaR=essaR+sr(ir)*FA
		 ettvR=ettvR-sr(ir)*HV*T13*ver
		 essvR=ettvR-sr(ir)*FV*T13*ver
	         ettmR=ettmR-sr(ir)*HV*t14*ver
	         essmR=essmR-sr(ir)*FV*t14*ver

		 vetar=vetar-sr(ir)*HV*w1
		 vesar=vesar-sr(ir)*FV*w1
		 getar=getar+sr(ir)*HV*dlr(ir)/dldop
		 gesar=gesar+sr(ir)*FV*dlr(ir)/dldop
102	      continue         !fin del do en componentes r Zeeman

	      ESAR=2.*ESAR
	      VESAR=2.*VESAR
	      GESAR=2.*GESAR
	      ESSAR=2.*ESSAR
	      ESSVR=2.*ESSVR
	      essmr=2.*essmr

	return
	end

