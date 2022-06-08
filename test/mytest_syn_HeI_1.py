import hazel
import matplotlib.pyplot as plt
import sys

mod = hazel.Model(working_mode='synthesis',atomf='helium.atom',verbose=0) 

mod.add_spectral({'Name': 'sp1', 'Wavelength': [10826, 10833, 150], 
	'topology': 'ch1','LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0]}) #what means 1 here??

#common arguments:
arg={'spectral region': 'sp1', 'atom':'helium','line': '10830', 
	'wavelength': [10826, 10833],'ref frame': 'LOS','coordB': 'spherical'}
#'reference frame': 'line-of-sight' or 'vertical'
#'coordB' : 'spherical' or 'cartesian' 

ch1=mod.add_chrom({'name': 'ch1','height': 3.0,**arg})
mod.setup() 

f, ax = plt.subplots(nrows=2, ncols=2)	;ax = ax.flatten()
for j in range(5):
	#PARS:(Bx or B,By or thB,Bz or phB,tau,v,deltav,beta,a).
	#OPT kwds:ff(default 1.0),j10(default np.zeros(4))	
	ch1.set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.01,0.01,0.01,0.01])
	mod.synthesize() 
	ax=mod.iplot_stokes(ax,'sp1')

plt.show()
