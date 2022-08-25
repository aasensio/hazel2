import hazel
import matplotlib.pyplot as pl
import numpy as np
import sys


nx=150
s0,sx=np.ones(nx),np.zeros(nx)
#----------------------------------------------------------------------------------
mo = hazel.Model(working_mode='synthesis',atomf='sodium_hfs.atom',verbose=0) 

cdic={'ref frame': 'LOS'}#common args
chs,txt=mo.add_Nchroms(['c0','c1'],cdic,hz=[0.,0.]) #return chs objects and tags
#chs,txt=mo.add_Nchroms(['c0','c1','c2'],cdic,hz=[0.,0.,0.]) #return chs objects and tags

s1=mo.add_spectrum('s1', atom='sodium',linehazel='5895',wavelength=[5892, 5899, nx], 
	topology='c0->c1',los=[0.0,0.0,90.0], boundary=[s0,sx,sx,sx])	;mo.setup() 

#----------------------------------------------------------------------------------

pl.close()	;f, ax = pl.subplots(nrows=2, ncols=2)	;ax = ax.flatten()
for j in [1,2,3]:
	#PARS:(Bx or B,By or thB,Bz or phB,tau,v,deltav,beta,a).j10=[t1,...tn], nbar=[t1,...tn], omega=x
	chs[0].set_pars([20.*j,10.*j,10.*j,2.,0.,8.,0.1,0.02],j10=[0.0,0.0],j20f=1.0)
	chs[1].set_pars([20.*j,10.*j,10.*j,2.,0.,8.,0.1,0.02],j10=0.0,j20f=[1.0,1.0])#
	#chs[2].set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.0,0.0])
	mo.synthesize() 
	ax=mo.iplot_stokes(ax,'s1')
pl.show()

f,ax2=mo.plot_coeffs('s1')

#ch1.set_pars([20.,10.,10.,3.,0.,7.,1.,0.02],j10=[0.,0.,0.1,0.1])
#ch2.set_pars([20.,80.,10.,3.,6.,7.,1.,0.2],j10=[0.,0.,0.3,0.2])
#mo.synthesize() 
#ax=mo.plot_stokes('s1')  #NON-interactive plot (all plotting things are inside) 

'''

atompol, magopt, stimem = 1, 1, 1  #
nocoh = 0  #=0 for including all coherences, set to level number to deactive cohs in it
dcol=np.asarray([0.0,0.0,0.0])  #D^(K=1 and 2)=delta_collision, D^(K=1),D^(K=2)
nbar = np.asarray([0.0,0.0,0.0,0.0])     #these four dimensions are maximum number of transitions
omega = np.asarray([0.0,0.0,0.0,0.0])   #not stokes components
#--------intensity profile to substitute constant I0[:,0] above----
I0 = np.zeros((nLambda,4))   #init boundary cond. CAREFUL:simplified slab uses Io=cte
#temp,dlamd=get_T(DoppWidth1,l0s[tran-1],atw=22.9897,vmic=vmic) #get temperature from thermal doppler width
#I0[:,0],ww = parametric_prof3(lamAxis-0-0,dep=0.9,uu=0.,amp=varlist[kk],alpha=dlamd,
#    gamma=0.2,rel=0.8,ratrel=0.99,sat=12,kind='2') #kind=2 requires dlamd
I0[:,0]=I0back 
#-------------------------
#synthesize
synMode = 5 #Opt.thin (0), slab no-MO (1), M-E (2), slab DELOPAR (3), simplified slab (4), exact slab (5)
'''