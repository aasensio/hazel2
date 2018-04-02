C CONHSRA construida con el programa CONTINUO
C calcula el continuo del HSRA por medio de un ajuste a un polinomio
c de tercer grado.El error maximo es de .125%

	real*4 function conhsra(x)

	implicit real*4 (a-h,o-z)
	real*4 c1(11),c2(11),c3(11),c4(11),c5(11),c6(11),c7(11)
	data c1/-4.906054765549e13,1.684734544039e11,1.507254517567e7
     &       ,-7561.242976546,0.,0.,0.,0.,0.,0.,0./
	data c2/-4.4650822755e14,6.1319780351059e11,-9.350928003805e7
     &       ,0.,0.,0.,0.,0.,0.,0.,0./
	data c3/-1.025961e15,1.3172859e12,-3.873465e8
     &       ,46486.541,-2.049,0.,0.,0.,0.,0.,0./
	data c4/4.861821e15,-2.2589885e12,4.3764376e8
     &       ,-39279.61444,1.34388,0.,0.,0.,0.,0.,0./
c	data c5/4.353082e17,-1.261996e14,5.31125e9,908890.899,43.24361
c     &       7,-2.0462285d-2,1.0018226758d-6,0.,0.,0.,0./
	data c5/1.758394e15,-3.293986e11,1.6782617e7,0.,0.,0.
     &       ,0.,0.,0.,0.,0./
	data c6/1.61455557e16,-6.544209e12,1.0159316e9
     &       ,-70695.58136,1.852022,0.,0.,0.,0.,0.,0./
	data c7/7.97805136e14,-1.16906597e11,5.315222e6
     &       ,-4.57327954,-3.473452d-3,0.,0.,0.,0.,0.,0./

	if(x.lt.3644.15)then
	    conhsra=c1(1)+x*(c1(2)+x*(c1(3)+x*(c1(4)+x*c1(5))))
	else if(x.lt.3750.)then
	    conhsra=c2(1)+x*(c2(2)+x*(c2(3)+x*(c2(4)+x*c2(5))))
	else if(x.lt.6250.)then
	    conhsra=c3(1)+x*(c3(2)+x*(c3(3)+x*(c3(4)+x*c3(5))))
	else if(x.lt.8300.)then
	    conhsra=c4(1)+x*(c4(2)+x*(c4(3)+x*(c4(4)+x*c4(5))))
	else if(x.lt.8850.)then
	    conhsra=c5(1)+x*(c5(2)+x*(c5(3)+x*(c5(4)+x*c5(5))))
	else if(x.lt.10000.)then
	    conhsra=c6(1)+x*(c6(2)+x*(c6(3)+x*(c6(4)+x*c6(5))))
	else
	    conhsra=c7(1)+x*(c7(2)+x*(c7(3)+x*(c7(4)+x*c7(5))))
	end if
	return
	end
