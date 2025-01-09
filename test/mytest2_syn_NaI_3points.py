import hazel
import matplotlib.pyplot as pl
import numpy as np
import sys

nx=150
p1,p0=np.ones(nx),np.zeros(nx) 

#----------------------------------------------------------------------------------
m1 = hazel.Model(mode='synthesis',atomfile='sodium_hfs.atom',apmosekc='1110')

cdic={'ref frame': 'LOS'}#common args to all chromospheres

chs,topo = m1.add_funcatmos(6,cdic,hzlims=[0.,1500.])#add functional atmosphere 
s1=m1.add_spectrum('s1', atom='sodium',linehazel='5895',wavelength=[5894, 5897, nx], 
	topology=topo,los=[0.,0.,90.],boundary=[p1,p0,p0,p0])	;m1.setup() 

#----------------------------------------------------------------------------------
dlims={'B1':[350,100], 'B2': [89,90], 'B3':[44,49],\
	'tau':[3.,0.05],'v':[0,4],'deltav':[4,7],'a':[0.2,0.1] ,\
	'j10':[0.01,0.02],'j20f':[1,1.5],'beta':[1,1]} #...'ff':[1,1],'nbar':[1,1]}

pkws={'plotit':9,'nps':3,'var':'mono','method':1}
hz=m1.set_funcatm(dlims,hztype='parab',orders=4,**pkws) #set atm pars following a given-order function

m1.synthesize(plot='s1') #frac=True
#m1.plot_coeffs('s1') #,coefs=['epsv','etai','etaq','etav'],scale=2)

#EXAMPLES MUTATES:
#mo,kk=m1.mutates('s1', apmosekc='1110',B1=[46.,48.],j20f=[1.1,1.5],pkws=pkws)
#dd={'B1':[1,34.],'j10'=[1,0.02],'method':'Emissivity'}
#m1.mutates('specname', atmpar1=[layernumber,value],j10=[layernumber,value],apmosekc='value',parsdic=dd)
#m1.mutates('s1', B1=[10,10.,1000],j10=[0,0.01,0.1])
#m1.mutates('s1', apmosekc='0110')
#m1.mutates('s1', apmosekc='0110',B1=[1,40.],j10=[0,0.02],bylayer=True)
#mo,kk=m1.mutates('s1', apmosekc='1111',B2=[1,0.],j10=[0,0.1],bylayer=True)
#mo,kk=m1.mutates('s1', apmosekc='1110',B3=[46.,48.],j20f=[1.1,1.5],pkws=pkws)
#mo2,kk=m1.mutates('s1', apmosekc='1110',v=[4.,8.],j10=[0.01,0.07],pkws=pkws)

#mo2.compare_experiments(mo,'s1')

#m1.fractional_polarization(s1)

#m1.reshow('all')#m1.reshow('1')

#ms=np.zeros(3, dtype=object)#gives a vector of pointers to python objects
#ms[0]=m1 try this to store all models of an experiment
