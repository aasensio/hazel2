c inicia_pefrompgt calcula la pe a partir de Pg y T usando Saha 
c solo para le H
    	subroutine inicia_pefrompgt(t,pg,pe)
 	real*4 t,pg,pe,nu,saha,aaa,bbb,ccc,ybh  	
    	
    	
    	nu=0.9091       ! only Hydrogen is ionized
	saha=-0.4771+2.5*alog10(t)-alog10(pg)-(13.6*5040./t)
	saha=10**saha
	aaa=1.d0+saha
	bbb=-(nu-1.)*saha
	ccc=-saha*nu

	ybh=(-bbb+sqrt(bbb*bbb-4.*aaa*ccc))/(2.*aaa)! ionization fraction
	
	pe=pg*ybh/(1.+ybh)
c       pe=pg*ybh

	return
	end