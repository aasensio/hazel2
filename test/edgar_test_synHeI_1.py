import hazel
import matplotlib.pyplot as plt
import sys

mo = hazel.Model(working_mode='synthesis',atomf='helium.atom',verbose=0) #'helium.atom' is default

ckeys={'ref frame': 'LOS'} # common args .  {'spectral region': 'sp1'} does not exist any more
#'ref frame': 'LOS' or 'vertical'(default)
#'coordB' : 'spherical'(DEFAULT) or 'cartesian' #keyword 'coordB' was 'coordinates for magnetic field vector'

#ch1=mo.add_chrom({'name': 'ch1','height': 1.0,**ckeys})  #OLD
#ch2=mo.add_chrom({'name': 'ch2','height': 2.0,**ckeys})

chs,tags=mo.add_Nchroms(['c1','c2'],ckeys,hz=[1.,3.]) #return a list with chs objects and their names

#the association between spectrum(or spectral region) with atmospheres in topology is now made here
s1=mo.add_spectrum('s1', atom='helium',linehazel='10830',wavelength=[10826, 10833, 150], 
	topology='c1->c2',los=[0.0,0.0,90.0], boundary=[1.0,0.0,0.0,0]) #i0fraction=1.0
'''
linesSIR='number' or '[number]' for adding lines to calcualte photospheres with SIR in the model
'atm_window': [10826, 10833] when omitted it takes range in 'Wavelength'
OLD: mo.add_spectral({'name':'sp1','atom':'helium','lineHazel': '10830','Wavelength': [10826, 10833, 150], 
	'topology':'c1->c2','LOS':[0.0,0.0,90.0],'Boundary condition': [1.0,0.0,0.0,0]}) 
here the 1 of boundary condition is the fraction of I0Allen that you assume for the background
boundary condition and for the normalization by i0 (1 means i0Allen at disk center). 


If a user background boundary condition is specified (via other keyword)...

Actually the normalization to the continuum intensity only requires the first number of the vector
boundary, while the other three zeros are always redundant, except when one wants to specify the
Stokes wavelength profile in the background. So:
1)we define an Icont normalization keyword containing the fraction of i0Allen to which we want to
normalize the output Stokes vector.
2)we adapt the boundary keyword to include wavelength-dependent profiles that can be defined by the user
in the present python program
3)we create a boundary subroutine to be called in synthesize of chromosphere.py where 
the above boundary conditions are translated to the Hazel boundaryIn parameter in chromosphere.py and to the
normalization to the I of the continuum. 

'''
mo.setup() 
'''
plt.close()
f, ax = plt.subplots(nrows=2, ncols=2)	;ax = ax.flatten()
for j in range(5):
	###PARS:(Bx/B,By/thB,Bz/phB,tau,v,deltav,beta,a).OPT kwds:ff(default 1.0),j10(default np.zeros(4))	
	chs[0].set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.,0.,0.2,0.])#,j10=[1e-2,2e-2,3e-2,4e-2])
	chs[1].set_pars([20.*j,10.*j,10.*j,3.,0.,8.,1.,0.02],j10=[0.,0.,0.,0.])#,j10=[1e-2,2e-2,3e-2,4e-2])
	mo.synthesize() 
	ax=mo.iplot_stokes(ax,'s1') #Interactive plot for loops
plt.show()

'''
chs[0].set_pars([20.,10.,10.,3.,0.,7.,1.,0.02],j10=[0.,0.,0.1,0.1])
chs[1].set_pars([20.,80.,10.,3.,6.,7.,1.,0.2],j10=[0.,0.,0.3,0.2])
mo.synthesize() 
ax=mo.plot_stokes('s1')  #NON-interactive plot (all plotting things are inside) 

f,ax2=mo.plot_coeffs('s1')

#TBD:

#pasa a sodio con dos capas y una iluminacion de frontera correcta.

#anade todos los cambios que hiciste en hazel1: 
#calculo de coeffs opticos,keywords, SEE, etc.
#testea calculo multiterm como en la GUI (archivo atomico?).
#reproduce resultados modelo.
#introduce rutinas fortran de transporte.
#Corrige EVOP.
#testea inversion : pasar j10 a componentes.
#crea keyword palette as made with "atmos window" in model.py.
#permite el calculo de multiterm y multiatom por separado para superponer resultado.
#añade teoria B.
#anade teoria C.

#prueba la gui.

#el "wavelength" de add_chromosphere podriamos cogerlo de add_spectral

#DONE:

#----------------COMMIT 1----------------------
#tras compilar los archivos fuente de hazel en directorio hazel se manda una copia a :
#/opt/anaconda3/envs/hazel_38/lib/python3.8/site-packages/hazel-2018.9.22-py3.8-macosx-10.9-x86_64.egg/hazel
#si intentas cambiar archivos en este ultimo directorio no sirve de nada porque se sobreescribibiran de nuevo
#al compilar desde el directorio hazel.

#En la linea 704 de model.py tienes "atom='sodium' " hardcoded y eso da error porque 
#el IF de la rutina add_active_line no pilla 10830 como una linea activa posible y salta en
# line_to_index[self.active_line] @chromosphere.py . 
#Asi que cambie la instancia de la clase Hazel_atmosphere:
#self.atmospheres[atm['name']] = Hazel_atmosphere(working_mode=self.working_mode, name=atm['name'], atom='sodium')
#Por esto mas generico:
#self.atmospheres[atm['name']] =Hazel_atmosphere(working_mode=self.working_mode, name=atm['name'], atom=atm['atom']) #EDGAR, quite atom='helium'

#En los archivos de ejemplo se te olvido una coma entre el comentario de Bz y tau y confunde.

#cambie el sistema de coordenadas para B por defecto a esfericas en model.py:
#self.atmospheres[atm['name']].coordinates_B = 'spherical'
#y le cambie el nombre al keyword 'coordinates for magnetic field vector'
#por 'coordB' tambien en model.py

#Introduje j10 como un vector de ntrans y todo lo necesario para que introducir
#j10 correctamente para cada atmosfera. 
#Luego habra que ver si la inversion se
#hace correctamente porque los parametros de la clase hazel atmopshre en 
#chromosphere.py no estan puestos como vector y ademas la inversion parece requerir
#numeros, no vectores.

#meti nueva rutina read_model_atom_ok en io_py.90 para leer archivos atomicos
# .atom genericos desde el directorio data .

#anadi alias y return mod.chromospheres dictionary a add_chromosphere, para 
#evitar hacer a1=mod.atmospheres['ch1'] o usar un nombre tan largo al 
#hacer set_parameters

#anadi alias a set_parameters como set_pars.

#quite el synthesisze despues de setup porque parece que no hace nada
#meti el synthesize despues de set_parameters dentro de set_parameters, y para
#ello le paso mod como argumento, pero esto es ineficiente y lo volvi a dejar fuera.
#hay 2 sintesize , uno para sintetizar el rtansporte en chromopshere.py 
#y otro para convolucionar y otras cosas en model.py.

#modifique rutina init() para aceptar verbose como parametro e iguale el verbose_mode
#de las rutinas core de hazel al verbose de las rutinas python.

#cambie el kwyword reference frame por uno mas corto y mejora la eficiencia
#de los IFs cuando se setea en add chromosphere

#anadi rutina de plot_stokes en model.py
# 	for i in range(4):ax[i].plot(mod.spectrum['sp1'].stokes[i,:])
#ahora es:
#	ax1=mod.plot_stokes(ax1,'sp1')
#----------------------- COMMIT 2

#corrected error in io_py.f90, where the module vars containing the atom strcuture was
#not being read.

#archivo atomico introducido dinamicamente(programaticamente) por python.Para esto hay que
#parchear el paso de cadena python a vector de caracterres en f90 a traves de wrapper pyx.
#el paso del nombre del archivo se hace a traves de la rutina init cuando se inicia el 
#el modulo Hazel en model.py. El path del archivo se presupone en hazel/data.

#...

#----------------------- COMMIT 5
#adicion de comentarios de desarrollador varios
#corregido un bug en la deteccion del tipo de "coordinates for magnetic field vector" en model.py

#-----------------------COMMIT 6
#Moved atom keyword from add_chromosphere to add_spectral. For this, we define a dictionary of dictionaries with 
#the atoms and lines and indexes of the spectral lines to choose among in model.py. The defintion of 
#line_to_index is then done reading from the dictionary from add_spectral, and in order to pass it to the 
#atmosphere object and to Hazel synthesize in chromosphere.py we do the following in the radiative transfer 
#logic of synthesize_spectral_region: self.atmospheres[atm].line_to_index=self.line_to_index.  

#-----------------------COMMIT 7
#ECR: moved line keyword and variables to add_spectral and to spectrum object from hazel and SIR atmospheres.

#-----------------------COMMIT 8
'''
ECR: Made necessary changes to change the order of add_spectral and add_chromosphere, 
finally disentangling atmopsheric aspects from spectral ones and reducing the verbosity 
and complexity of the code. To this aim, add_active_line is now done for every atmosphere inside 
add_spectral(right after adding all atmospheres of the topology).Before, both add_active_line and 
the block before appeared repeated for every kind of adding atmosphere routine,preventing to change
 the order of add_spectral and add_chromosphere. Now n_chromospheres can be known and calculated 
 from add_spectral, before setup, which is the requirement for starting to program the storage of 
 all optical coefficients.
'''

#-----------------------COMMIT 9
'''
ECR: Created add_spectrum in model.py to substitute add_spectral. Now it works with keywords instead of with input dictionaries.
Remove "name" from input dictionaries of add_spectrum , leaving it as fixed parameter.
Add additional reduced keyword 'boundary' instead of 'boundary condition' in add_spectrum.
Create containers for optical coefficients in spectrum object.
Modify fortran routines to extract optical coefficients of every chromosphere to Python. 
Check transfer in synthesize_spectral_region.
Add routine for adding N chromospheres in one line.
Add plotting routines for optical coefficients.


#-----------------------COMMIT 10
Fix bug when creating a model with other atom different than helium. Now there is a dictionary with 
multiplets at the beginning of model.py, such that all atom-related information is together and more clear.
Change to sodium atom creating test program.
Improved efficiency and readibility of magnetic field transformations in synthesize of chromosphere.py calling Hazel.
Addition of subroutines cartesian_to_spherical, spherical_to_cartesian, and los_to_vertical.
Addition of get_B_Hazel and simplification of the treatment of the magnetic field reference and coordinates
in add_parameters and synthesize inside chromosphere.py.
CHeck that j10 is correctly introduced(OK).
Check that effectively linear anisotropy is introducing J10 in the SEE of sodium(OK).
Reduced parsInput to parsIn names in synthesize at chromosphere.py. 
Minor simplifications in synthesize_spectral_region.


#-----------------------COMMIT 11
Definition of the dictionary atms_in_spectrum in model.py to link the name of the spectra to 
the topologies, such that the outer loop and the if clause in synthesize_spectral_region can be
avoided. Now the loop does not need to search over all atmospheres but just go through the ones
belonging to that spectral region of interest.
Fix bug when storing optical coefficients in topologies of the kind c0->c1+c2. Solved 
using index n+k in synthesize_spectral_region. Simplification of the routine synthesize_spectral_region
into a new routine called synthesize_spectrum, where unnecessary operations and code lines are removed.
Detected bug in continuum of Stokes I (being 2.0 instead of 1.0) in topologies of the kind c0->c1+c2. 
 This actually happens when the user forget to specify the ff's for the atmospheres
in set_pars or set_parameters. Possible bug in old routine normalize_ff, it seems not working for synthesis.
To fix these things we add check_filling_factors routine, which assures that the sum of ff's in 
every layer with subshells amounts to 1. If not, it assumes isocontribution: ff=1/natms with nsub the number 
of sub-atmospheres in the layer (e.g.,if nsub=2, ff=0.5). 


#-----------------------COMMIT 12(pending) : just before SUMMER holidays:
New procedure for introducing boundary conditions with spectral dependence. 
Introduction of i0fraction keyword as the fraction of the true 
background continuum intensity to the Allen continuum.
i0fraction is included just in case we want to consider it as inversion parameter in the future, but 
 in general it shall be avoided (is set to 1.0 by default). When defining boundary as the ratio between
 the true physical background Stokes profiles and the I0Allen, i0fraction is already the first
 spectral value of the boundary intensity condition, being 1 to represent i0Allen or lower for different
 continuum backgrounds. 
Now the boundary keyword accept either 4 scalars that will be broadcasted as constants in 
wavelength for every Stokes, or directly 4 spectral profiles. 
These changes shall allow to control the boundary spectral dependence and the relative level of
continuum intensity that shall be used to normalize, allowing treating solar
locations where I0Allen does not fit the profiles. In addition to the changes in model.py, 
spectrum.py, and chromospehre.py to simplify and generalize the set up of the boundary condition, I 
also added the check_key routine to simplify the setup of default keyword values 
in the (deprecated) add_spectral subroutine.

I changed the number of dimensions of j10 , nbar, and omega to ntrans, thus avoiding a hardcoded
 value (it was hardcoded to 4). The flow is :read fortran atom%ntrans from atom file (from io_f90.py) -> 
 init() routine in hazel_py.f90 -> init() in hazel_code.pyx -> call hazel_code._init in model.py, 
 where self.ntrans is set and later passed to set_chromosphere (from main python program or from
 file with use_configuration), which passes it via input parameter to initialization of Hazel chromosphere
object, where dimensions of j10,omega and nbar are set to ntrans.  


CAUTION: in helium.atom el valor del nbar reduction factor de la transicion 2 esta puesto 
a 0.2 en vez de a 1.0. Verificar si esto se sobreescribe luego en el codigo y preguntar a Andres
 por qué debe estar asi.
--



Reintroduce  UseAtomicPol (which is set always to 1 in hazel_pz.f90) and similar parameters.
Start to add all your Hazel1 changes to Hazel2.
clean commments

'''





