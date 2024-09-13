import hazel
import matplotlib.pyplot as pl
import numpy as np
import sys

nx=150
p1,p0=np.ones(nx),np.zeros(nx) 

#----------------------------------------------------------------------------------
m1 = hazel.Model(mode='synthesis',atomfile='sodium_hfs.atom',apmosekc='1110')

cdic={'ref frame': 'LOS'}#common args to all chromospheres

chs,topo = m1.add_funcatmos(6,cdic,hzlims=[0.,1500.0])#add functional atmosphere 
s1=m1.add_spectrum('s1', atom='sodium',linehazel='5895',wavelength=[5894, 5897, nx], 
	topology=topo,los=[0.,0.,90.],boundary=[p1,p0,p0,p0])	;m1.setup() 

#----------------------------------------------------------------------------------
dlims={'B1':[350,100], 'B2': [89,90], 'B3':[44,49],\
	'tau':[3.,0.05],'v':[0,4],'deltav':[4,7],'a':[0.2,0.1] ,\
	'j10':[0.01,0.02],'j20f':[1,1.5],'beta':[1,1]} #...'ff':[1,1],'nbar':[1,1]}

pkws={'plotit':9,'nps':3,'var':'mono','method':1}
hz,chs,dat=m1.set_funcatm(dlims,hztype='parab',orders=4,**pkws) #set atm pars following a given-order function

COMMIT:
Creation of subroutines to synthesize optically-thick lines better:
1) add_funcatmos(), set_funcatmos() and PolyFx() to easily and succintely add and set parametric 
 optically-thick atmosphere made of N Hazel atmospheres with physical quantities
varying as different analytical functions (polynomials,exponential,random,...). 
2) plot_funcatmos() and plot_PolyFx() to plot the variations with height of all atmospheric
parameters and compare them with polynomials of different order.
3) fix_point_polyfit_fx() to create a polynomial fitting certain points but
passing exactly through some control points (typically the boundaries).
4) fun_minT() to easily create parametric variations of deltav mimicking a minimum 
of temperature.
5)get_Tpars() to obtain transformations between temperature, doppler broadening, and deltav.
Adding here a new dictionary of atomic weights in model. 
6) check_B_vals() to check that the magnetic field values are in the correct limits
for every coordinate system.Add new dictionary of limiting values for each physical quantity.
7)Allow the height axis to be non-linear to check the effects of sampling some
 optical depths more than others.



#m1.synthesize(plot='s1') #frac=True
#m1.plot_coeffs('s1')#,coefs=['epsv','etai','etaq','etav']),scale=2

#dd={'B1':[1,34.],'j10'=[1,0.02],'method':'Emissivity'}
#m1.mutates('specname', atmpar1=[layernumber,value],j10=[layernumber,value],apmosekc='value',parsdic=dd)
#m1.mutates('s1', B1=[10,10.,1000],j10=[0,0.01,0.1])
#m1.mutates('s1', apmosekc='0110')
#m1.mutates('s1', apmosekc='0110',B1=[1,40.],j10=[0,0.02])

#m1.fractional_polarization(s1)

