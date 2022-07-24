import hazel
import matplotlib.pyplot as plt
import sys

mod = hazel.Model(working_mode='synthesis',atomf='helium.atom',verbose=0) 

arg={'ref frame': 'LOS'}#common arguments
ch1=mod.add_chrom({'name': 'ch1','height': 1.0,**arg})
ch2=mod.add_chrom({'name': 'ch2','height': 3.0,**arg})

mod.add_spectral({'Name': 'sp1', 'atom':'helium','lineHazel': '10830', 'Wavelength': [10826, 10833, 150], 
	'topology': 'ch1->ch2','LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0]})

mod.setup() 

plt.close()	;f, ax = plt.subplots(nrows=2, ncols=2)	;ax = ax.flatten()
for j in range(5):
	#PARS:(Bx or B,By or thB,Bz or phB,tau,v,deltav,beta,a).
	#OPT kwds:ff(default 1.0),j10(default np.zeros(4))	
	ch1.set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.01,0.01,0.01,0.01])
	ch2.set_pars([20.*j,10.*j,10.*j,3.,3.,8.,1.,0.02],j10=[0.01,0.01,0.01,0.01])
	mod.synthesize() 
	ax=mod.iplot_stokes(ax,'sp1')

plt.show()


#ch1.set_pars([20.,10.,10.,3.,0.,7.,1.,0.02],j10=[0.,0.,0.1,0.1])
#ch2.set_pars([20.,80.,10.,3.,6.,7.,1.,0.2],j10=[0.,0.,0.3,0.2])
#mod.synthesize() 
#ax=mod.plot_stokes('sp1')  #NON-interactive plot (all plotting things are inside) 
