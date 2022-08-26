import hazel
import matplotlib.pyplot as pl
import numpy as np
import sys

nx=150
s0,sx=np.ones(nx),np.zeros(nx) 
#----------------------------------------------------------------------------------
mo = hazel.Model(mode='synthesis',atomfile='sodium_hfs.atom', apmosenc='1110',dcol=[0.,0.,0.])

cdic={'ref frame': 'LOS'}#common args to all chromospheres
chs,txt=mo.add_Nchroms(['c0','c1'],cdic,hz=[0.,0.]) #return chs objects and tags

s1=mo.add_spectrum('s1', atom='sodium',linehazel='5895',wavelength=[5892, 5899, nx], 
	topology='c0->c1',los=[0.,0.,90.],boundary=[s0,sx,sx,sx])	;mo.setup() 

#----------------------------------------------------------------------------------

pl.close()	;f, ax = pl.subplots(2,2)	;ax = ax.flatten()
for j in [1,2,3]:
	#PARS:(Bx or B,By or thB,Bz or phB,tau,v,deltav,beta,a).j10=[t1,...tn],j20f(omega)=[v1,...vn],nbar=[x1,...xn]
	chs[0].set_pars([20.*j,10.*j,10.*j,2.,0.,8.,0.1,0.02],j10=0.0,j20f=1.0)#,nbar=1.0) 
	chs[1].set_pars([20.*j,10.*j,10.*j,2.,0.,8.,0.1,0.02],j10=0.0,j20f=1.0)
	
	mo.synthesize(method='Emissivity')#'Emissivity' method='EvolOp',muAllen=1.0
	ax=mo.iplot_stokes(ax,'s1')

pl.show()
f2,ax2=mo.plot_coeffs('s1')

'''
newmo=mo.mutation({},apmosenc='1010',compare=True)

develop exp_mutation: le pasas el modulo del experimento anterior, el s1, y las cromosferas,
y luego le dices que cambie unos pocos parametros que tu quieras y compre el resultado con el anterior.

'''