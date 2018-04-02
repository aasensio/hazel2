	subroutine MVOIGTC(a,ver,t13,t14,wc,dldop,H,vetar,
     &      HA,ettvr,ettmr)

	REAL*4 WC
	common/piis/piis
		
	call voigt2(a,ver,H,F)
	HV=-2.*ver*H+4.*a*F           !parcial de H con respecto a v_doppler
	HA=-2.*(piis-a*H-2.*ver*F)    !derivada con a de H
	ettvR=-HV*t13*ver             !derivada de H con resp. a T por 
				      !cambios en anchura doppler
	ettmR=-HV*t14*ver	      !derivada de H con mic
	vetar=-HV*wc/dldop

	return
	end

