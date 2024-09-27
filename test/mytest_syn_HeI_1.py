import hazel
import matplotlib.pyplot as plt
import sys

#----------------------------------------------------------------------------------
mo = hazel.Model(working_mode='synthesis',atomf='helium.atom',verbose=0) 

cdic={'ref frame': 'LOS'}#common args
chs,txt=mo.add_Nchroms(['c1','c2'],cdic,hz=[1.,3.]) #return chs objects and tags

s1=mo.add_spectrum('s1', atom='helium',linehazel='10830',wavelength=[10826, 10833, 150], 
	topology='c1->c2',los=[0.0,0.0,90.0], boundary=[1.0,0.0,0.0,0])	;mo.setup() 

#----------------------------------------------------------------------------------

plt.close()	;f, ax = plt.subplots(nrows=2, ncols=2)	;ax = ax.flatten()
for j in range(5):
	#PARS:(Bx or B,By or thB,Bz or phB,tau,v,deltav,beta,a).
	chs[0].set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.01,0.01,0.01,0.01])
	chs[1].set_pars([20.*j,10.*j,10.*j,3.,2.,8.,1.,0.02],j10=[0.01,0.01,0.01,0.01])
	mo.synthesize() 
	ax=mo.iplot_stokes(ax,'s1')
plt.show()


f,ax2=mo.plot_coeffs('s1')

#ch1.set_pars([20.,10.,10.,3.,0.,7.,1.,0.02],j10=[0.,0.,0.1,0.1])
#ch2.set_pars([20.,80.,10.,3.,6.,7.,1.,0.2],j10=[0.,0.,0.3,0.2])
#mo.synthesize() 
#ax=mo.plot_stokes('s1')  #NON-interactive plot (all plotting things are inside) 
