import hazel
import matplotlib.pyplot as pl
import numpy as np
import sys


nx=150
s0,sx=np.ones(nx),np.zeros(nx) 
'''
edic={'Atompol':1,'MO effects':1,'Stim. emission':1, 'Kill coherences':0,'dcol':[0.,0.,0.]}
dcol=[0.,0.,0.],extrapars=edic,'helium.atom' and verbose =0 are default

ap-mo-se-nc --> atompol, magopt, stimem, nocoh = 1, 1, 1, 0  
#nocoh =0 includes all coherences, set to level number to deactive cohs in it

dcol=[0.0,0.0,0.0]  #D^(K=1 and 2)=delta_collision, D^(K=1),D^(K=2)

synMode = 5 #Opt.thin (0), DELOPAR (3),  EvolOp (5)
'''
#----------------------------------------------------------------------------------
m1 = hazel.Model(mode='synthesis',atomfile='sodium_hfs.atom',apmosekc='1110')

cdic={'ref frame': 'LOS'}#common args to all chromospheres
chs,txt=m1.add_Nchroms(['c0','c1'],cdic,hz=[0.,0.]) #return chs objects and tags

s1=m1.add_spectrum('s1', atom='sodium',linehazel='5895',wavelength=[5894, 5897, nx], 
	topology='c0->c1',los=[0.,0.,90.],boundary=[s0,sx,sx,sx])	;m1.setup() 

#----------------------------------------------------------------------------------

for j in [3,2,1]:#(Bx\B,By\thB,Bz\phB,tau,v,deltav,beta,a),ff=x,j10=[t1,...tn],j20f=[v1,...vn],nbar=[x1,...xn]
	chs[0].set_pars([20.*j,10.*j,10.*j,2.,0.,8.,0.1,0.02],j10=0.0,j20f=1.0)#,nbar=1.0) 
	chs[1].set_pars([20.*j,10.*j,10.*j,2.,3.,9.,0.1,0.02],j10=0.0,j20f=1.0)
	m1.synthesize(plot='s1') #frac=True

m1.plot_coeffs('s1')#,coefs=['epsv','etai','etaq','etav']),scale=2

#dd={'B1':[1,34.],'j10'=[1,0.02],'method':'Emissivity'}
#m1.mutates('specname', atmpar1=[layernumber,value],j10=[layernumber,value],apmosekc='value',parsdic=dd)
m1.mutates('s1', B1=[1,34.],j10=[0,0.01])
m1.mutates('s1', apmosekc='0110')
m1.mutates('s1', apmosekc='0110',B1=[1,40.],j10=[0,0.02])


#m1.fractional_polarization(s1)

'''
CHECK:
haz commit, duplica el fork, actualiza remotos, cambia documentacion poniendo ejemplos, y/o add 
self-explanatory help functions to every module and subroutine.
1) shapes of rho versus eta coeffs (is atompol affecting only etas?)
2) atompol is being acting in coeffs when j10=0???
3) units and amplitudes of coeffs for different LOS's, why q y u 4 orders of magnitude smaller than v? 
4)memory usage...may be you should add some garbage collector routines.
5)Add functional variations of many atmopshres and plot_coeffs_2D for many atmospheres









'''