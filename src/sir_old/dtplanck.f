c	dtplanck es una funcion que calcula la derivada de la funcion
c	de planck con respecto a la temperatura

	real*4 function dtplanck(t,lambda)
	real*4 t,lambda,c1,c2,d4,b,ex,dplnck

	c1=1.1910627e-5		!1./(2*h*c**2) con lambda en cm
	c2=1.43879
	d4=lambda*lambda*lambda*lambda
	b=dplnck(t,lambda)
	ex=c2/(lambda*t)
	dtplanck=(c2*b*b*d4*exp(ex))/(c1*t*t)

	return
	end
