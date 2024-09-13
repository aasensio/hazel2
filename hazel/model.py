from hazel.chromosphere import Hazel_atmosphere
from hazel.photosphere import SIR_atmosphere
from hazel.parametric import Parametric_atmosphere
from hazel.stray import Straylight_atmosphere
from hazel.configuration import Configuration
from hazel.io import Generic_output_file
from collections import OrderedDict
from hazel.codes import hazel_code, sir_code
from hazel.spectrum import Spectrum
from hazel.transforms import transformed_to_physical, physical_to_transformed, jacobian_transformation
import hazel.util
import numpy as np
import copy
import os
from pathlib import Path
import scipy.stats
import scipy.special
import scipy.signal
import scipy.linalg
import scipy.optimize
import warnings
import logging
import sys
import matplotlib.pyplot as plt #EDGAR: Im placing plotting routines here, but is a bit ugly


labdic = {'z1':r'$\mathrm{z \, [Mm]}$',
        'tt':r'$\mathrm{T\,[kK]}$','tit':r'$\mathrm{Temperature}$',
        'vdop':r'$\mathrm{V^{dop}_z \,[D.u.]}$',
        'b':r'$\mathrm{|B| \,[G]}$','tb':r'$\mathrm{\theta_B \,[Degrees]}$',
        'cb':r'$\mathrm{\chi_B \,[Degrees]}$',
        'frokq':r'$\{\rho}^K_Q/\rho^0_0$',
        'xx0':r'$\mathrm{\lambda-\lambda_0[{\AA}]}$',  
        'xx':r'$\mathrm{\lambda[{\AA}]}$',  
        'iic':r'$\mathrm{I/I_c}$','qi':r'$\mathrm{Q/I}$','ui':r'$\mathrm{U/I}$', 
        'vi':r'$\mathrm{V/I}$','qic':r'$\mathrm{Q/I_c}$','uic':r'$\mathrm{U/I_c}$', 'vic':r'$\mathrm{V/I_c}$', 
        'epsi':r'$\mathrm{\epsilon_I}$','epsq':r'$\mathrm{\epsilon_Q}$','epsu':r'$\mathrm{\epsilon_U}$',
        'epsv':r'$\mathrm{\epsilon_V}$','etai':r'$\mathrm{\eta_I}$','etaq':r'$\mathrm{\eta_Q}$',
        'etau':r'$\mathrm{\eta_U}$','etav':r'$\mathrm{\eta_V}$','rhoq':r'$\mathrm{\rho_Q}$',
        'rhou':r'$\mathrm{\rho_U}$','rhov':r'$\mathrm{\rho_V}$'
        }


def mylab(lab):
    return labdic.get(lab,lab) #return the input keyword string if lab is not in labdic

'''
def latex1(str):
    return r'$\mathrm{z \,{0} [Mm]}$'.format(str) 
'''
def exact_parabols(xax,P_ini,P_end,aval=0.2):
    #P_ini-->[p1x,p1y]
    #P_end-->[p2x,p2y]
    #define the function, the a parameter allows different parabols
    fn= lambda x,a,pa,pb : (a*(x-pb[0])+(pb[1]-pa[1])/(pb[0]-pa[0]))*(x-pa[0])+pa[1]
    return fn(xax,aval,P_ini,P_end)

def exp_3points(xax,p0,p1,p2):#needs 3 points to be fit 
    fn = lambda x,a,b,c : a + b*np.exp(c * x)
    return fn(xax,p0,p1,p2)


def exp_2points(xax, xl, yl):#fit with 2 points
    c=np.log(yl[0]/yl[1])/(xl[0]-xl[1])
    b=yl[1]*np.exp(-c*xl[1])
    return b*np.exp(c * xax)


__all__ = ['Model']

class Model(object):
    def __init__(self, config=None, mode='synthesis', atomfile='helium.atom',apmosekc='1110', dcol=None,
        extrapars=None, verbose=0, debug=False,rank=0, randomization=None, plotit=False, root=''):
        '''
        edic={'Atompol':1,'MO effects':1,'Stim. emission':1, 'Kill coherences':0,'dcol':[0.,0.,0.]}
        dcol=[0.,0.,0.],extrapars=edic,'helium.atom' and verbose =0 are default

        ap-mo-se-nc --> atompol, magopt, stimem, nocoh = 1, 1, 1, 0  
        #nocoh =0 includes all coherences, set to level number to deactive cohs in it

        dcol=[0.0,0.0,0.0]  #D^(K=1 and 2)=delta_collision, D^(K=1),D^(K=2)

        synMode = 5 #Opt.thin (0), DELOPAR (3),  EvolOp (5)
        '''

        #check/init mutable keywords at starting to avoid having memory of them among function/class calls
        if extrapars is None:extrapars={}
        if dcol is None:dcol=[0.,0.,0.]

        np.random.seed(123)
        if (rank != 0):
            return

        #EDGAR: dictionary of dictionaries with possible atoms and lines with their indexes for HAZEL atmospheres and spectra        
        #A variation of this dict to add more atoms and lines required to do similar changes in 
        #multiplets dict in general_atmosphere object (atmosphere.py). 
        #Now everything is encapsulated here in atomsdic.
        #From here we could also extract the number of transitions as len(self.multipletsdic[atom]).
        #However that would be a redundancy with respect to the data in the atom files,
        #and furthermore ntrans must be read at the very beginning before initializing j10, 
        #so the choice is to read ntrans from atom file: 
        #io_py.f90 ->hazel_py.f90 -> hazel_code.pyx -> init routine here in model.py

        self.atomfile=atomfile #memory for mutations

        self.atomsdic={'helium':{'10830': 1, '3888': 2, '7065': 3,'5876': 4},
            'sodium':{'5895': 1, '5889': 2}} 
        
        self.multipletsdic={'helium':{'10830': 10829.0911, '3888': 3888.6046, '7065': 7065.7085, '5876': 5875.9663},
                            'sodium':{'5895': 5895.924, '5889': 5889.95}}

        self.atwsdic={'helium':4.,'sodium':22.9897,'calcium':40.08} 


        self.apmosekc=apmosekc #memory for mutations, must be before get_apmosekcl
        self.dcol=dcol #memory for mutations, must be before get_apmosekcl

        #apmosekcL is List whose last element is other list with dcol
        self.apmosekcl=self.get_apmosekcl(apmosekc,dcol,extrapars) 


        #Dictio of minimum, default, and maximum values for all possible pars in Hazel atmosphere 
        #this could be conflicting with the use of "ranges", but such ranges seem to be applied only
        #in inversion and I dont see their actual set up to numeric meaningful values anywhere.
        #although 3 last pars are lists of ntr elements themselves, the limit is common to all of elements
        #use dmm in set_B_vals or make it global:  
        self.limB=4000.
        self.dmm={'Bx': [0.,100.,self.limB], 'By': [0.,100.,self.limB], 'Bz': [0.,100.,self.limB], \
            'B': [0.,100.,self.limB], 'thB': [0.,0.,180.], 'phB': [-360.,0.,360.], \
                'tau':[0.,1.,20.],'v':[-50.,0., 50.],'deltav':[0.5,2.,15.], \
                'beta':[0.,1.,10.],'a':[0.,0.1,10.],'ff':[0.,1.,1.], \
                'j10':[0.,0.,1.],'j20f':[0.,1.,1000.],'nbar':[0.,1.,1.]}
        


        self.plotit=plotit
        self.plotscale=3
        self.labelf1='1'
        self.labelf2='2'
        self.labelf3='3'
        self.f1=None
        self.f3=None
        self.ax1= None
        self.ax3=None
        self.ax5=None

        self.lock_fractional=None

        #synthesis methods to be implemented
        self.methods_dicT={0:'Emissivity',1:'Delo1',2:'Delo2',3:'Hermite',4:'Bezier',5:'EvolOp',6:'Guau'} 
        self.methods_dicS={'Emissivity':0,'Delo1':1,'Delo2':2,'Hermite':3,'Bezier':4,'EvolOp':5,'Guau':6} 
        self.methods_list=[ss for ss,tt in self.methods_dicS.items()] #list with only the names
        
        self.synmethod=5 #5 is default and can be changed by add_spectrum and /or by synthesize.
        #------------------------------------
        self.photospheres = []
        self.chromospheres = []
        self.chromospheres_order = []
        self.atmospheres = {}
        self.order_atmospheres = []        
        self.straylight = []
        self.parametric = []
        self.configuration = None
        self.n_cycles = 1
        self.spectrum = []
        self.spectrum = {}  #EDGAR:self.spectrum was initialized as [] and now as {}
        self.topologies = {}#[] EDGAR: before it was a list, now it is a dictionary
        self.atms_in_spectrum={} #EDGAR: of the kind -> {'sp1': string_with_atmosphere_order_for_sp1}
        
        #default mu where Allen continuum shall be taken for normalizing Stokes output
        #the actual value is set when calling synthesize_spectrum
        self.muAllen=1.0 

        self.straylights = []
        self.working_mode = mode
        self.pixel = 0
        self.debug = debug
        self.use_analytical_RF_if_possible = False
        self.nlte_available = False
        self.use_nlte = False
        self.root = root

        self.epsilon = 1e-2
        self.svd_tolerance = 1e-8
        self.step_limiter_inversion = 1.0
        self.backtracking = 'brent'
        
        self.verbose = verbose
        
        self.logger = logging.getLogger("model")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []

        ch = logging.StreamHandler()
        # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

        # Set randomization
        if (randomization is None):
            self.n_randomization = 1        
        else:
            self.n_randomization = randomization

        if (self.verbose >= 1):
            self.logger.info('Hazel2 Experimental')
        
        if ('torch' in sys.modules and 'torch_geometric' in sys.modules):
            if (self.verbose >= 1):
                self.logger.info('PyTorch and PyTorch Geometric found. NLTE for Ca II is available')
            self.nlte_available = True

        #We initialize pyhazel (and set self.ntrans) before calling add_chromosphere in use_configuration 
        #in order to setup self.ntrans before defining j10 length.
        #before these changes, the Hazel init was done after the next if..else.  
        self.ntrans=hazel_code._init(atomfile,verbose)   #EDGAR
        
        if (config is not None):
            if (self.verbose >= 1):
                self.logger.info('Using configuration from file : {0}'.format(config))
            self.configuration = Configuration(config)

            #EDGAR:n_chromospheres is set here.
            self.use_configuration(self.configuration.config_dict) 


    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d:
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d:
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)

    def __str__(self):
        tmp = ''
        for l, par in self.__dict__.items():
            if (l != 'LINES'):
                tmp += '{0}: {1}\n'.format(l, par)
        return tmp

    def get_apmosekcl(self,apmosekc,dcol,extrapars):
        '''
        Parameters in extrapars overwrite those in apmosekc.
        apmosekc and extrapars are mostly redudant on purpose. 
        apmosekc is much more compact but extrapars is there for who requires more readibility
        extrapars={'Atompol':1,'MO effects':1,'Stim. emission':1, 'Kill coherences':0,'dcol':[0.,0.,0.]}
        '''
        apmosekcl = [int(x) for x in apmosekc] #now is a list of atompol,magopt,stimem,nocoh
        for elem in apmosekcl[0:3]:
            if (elem != 0) and (elem!=1):raise Exception("ERROR: apmosekc first values can only be zeros or ones")
        #EDGAR: we should also check that the nocoh value does not go beyond number of levels(pending)

        for kk,keyw in enumerate(['Atompol','MO effects','Stim. emission','Kill coherences']):
            if (keyw in extrapars) and (extrapars[keyw]!=apmosekcl[kk]):
                apmosekcl[kk]=extrapars[keyw]        
                if (keyw!='Kill coherences') and (extrapars[keyw] != 0) and (extrapars[keyw]!=1):
                    raise Exception("ERROR: firsts parameters in extrapars can only be zeros or ones")    
        
        if ('dcol' in extrapars):apmosekcl.append(extrapars['dcol'])
        else:apmosekcl.append(dcol)

        return apmosekcl

    '''
    #EDGAR: this routine will not work if the backend is not compatible. 
    #By default I use MacOs backend, whose canvas does not have the window attribute.

    def move_figure(self,f, x, y):
        """Move figure's upper left corner to pixel (x, y)"""
        backend = matplotlib.get_backend()
        print(backend)
        if backend == 'TkAgg':
            f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
        elif backend == 'WXAgg':
            f.canvas.manager.window.SetPosition((x, y))
        else:
            # This works for QT and GTK
            # You can also use window.setGeometry
            f.canvas.manager.window.move(x, y)
    '''

    def setup_set_figure(self,scale=3):
        pscale=self.plotscale
        if scale!=pscale:pscale=scale

        plt.close('all')  
        plt.rcParams.update({'font.size': pscale*2.4})
        
        #fig, axs = plt.subplots(2, 2,figsize=(pscale*2,pscale*2))  
        self.f1, ax = plt.subplots(2, 2,figsize=(pscale*2,pscale*2),label=self.labelf1)  
        self.ax1 = ax.flatten()

        #start_x, start_y, dx, dy = self.f1.figbbox.extents
        #self.move_figure(self.f1,505,500)
        

    def toggle_visible(self, fig,visible=True):
        fig.set_visible(not fig.get_visible()) #if invisible, make it visible, or viceversa
        plt.draw()

    def remove_fig(self,fig_id):
        '''
        plt.close(X) X can be:
        *None*: the current figure
        .Figure: the given .Figure instance
        int: a figure number
        str: a figure name
        'all': all figures
        '''
        if (self.verbose >= 1):self.logger.info('Before removing_fig():',plt.get_figlabels()) #plt.get_fignums()
        
        plt.close(fig_id)
        if self.f1 is not None:
            self.f1.canvas.manager.close()
            self.f1= None
        
        if (self.verbose >= 1):self.logger.info('After removing_fig():',plt.get_figlabels()) #plt.get_fignums()
        
        plt.draw() #EDGAR:should not be removed??


    def fractional_polarization(self,sp,scale=3,tf=4,ax=None,lab=['iic','qi','ui','vi']): 
        '''
        Compares fractional and not fractional polarization.
        This routine is valid when fractional polarization is not implemented in synthesize_spectrum
        '''
        pscale=self.plotscale
        if scale!=pscale:pscale=scale

        plt.rcParams.update({'font.size': pscale*tf})
        if type(sp) is not hazel.spectrum.Spectrum:sp=self.spectrum[sp]

        if ax is None:
            f, ax = plt.subplots(nrows=2, ncols=2,figsize=(pscale*2,pscale*2))  
            ax = ax.flatten()

        labs=[r'$P/I_c$',r'$P/I$']
        
        for i in range(4):
            if i==0:
                ax[i].plot(sp.wavelength_axis, sp.stokes[i,:],color='r')
            else:
                l1,=ax[i].plot(sp.wavelength_axis, sp.stokes[i,:],color='r')
                l2,=ax[i].plot(sp.wavelength_axis, sp.stokes[i,:]/sp.stokes[0,:],color='k')
                
            ax[i].set_title(mylab(lab[i]))
            if i>1:ax[i].set_xlabel(mylab('xx'))
        
        #f.legend(tuple(lines), tuple(labs), loc=(0.1,0.1), bbox_to_anchor=(0.1, 0.3))
        plt.gcf().legend(tuple([l1,l2]), tuple(labs))#loc=(0.1,0.1), bbox_to_anchor=(0.1, 0.3))

        plt.tight_layout()
        plt.show()
        
        return 


    def plot_stokes(self,sp,scale=3,tf=4,fractional=False,lab=None): 
        '''
        Routine called by synthesize to plot stokes profiles either fractional 
        or normalized to continuum 
        '''
        if fractional:lab=['iic','qi','ui','vi']
        else:lab=['iic','qic','uic','vic']

        pscale=self.plotscale
        if scale!=pscale:pscale=scale

        plt.rcParams.update({'font.size': pscale*tf})
        if type(sp) is not hazel.spectrum.Spectrum:sp=self.spectrum[sp]

        if self.ax1 is None:setup_set_figure()

        for i in range(4): 
            if i ==0:
                self.ax1[i].plot(sp.wavelength_axis, sp.stokes[i,:])
            else:
                if fractional:self.ax1[i].plot(sp.wavelength_axis, sp.stokes[i,:]/sp.stokes[0,:])
                else:self.ax1[i].plot(sp.wavelength_axis, sp.stokes[i,:])
            
            self.ax1[i].set_title(mylab(lab[i]))#,size=8 + 0.7*pscale)
            if i>1:self.ax1[i].set_xlabel(mylab('xx'))#,size=8 +0.7*pscale)#,labelpad=lp)
        
        plt.tight_layout()
        plt.show()
        
        return 


    def plot_stokes_direct(self,sp,scale=3,tf=4,lab=None): 
        '''
        This plots stokes profiles normalized to continuum directly from spectrum object
        without considering fractional polarization 
        '''
        if lab is None:lab=['iic','qic','uic','vic']

        pscale=self.plotscale
        if scale!=pscale:pscale=scale

        plt.rcParams.update({'font.size': pscale*tf})
        if type(sp) is not hazel.spectrum.Spectrum:sp=self.spectrum[sp]

        if self.ax1 is None:setup_set_figure()

        for i in range(4): 
            self.ax1[i].plot(sp.wavelength_axis, sp.stokes[i,:])
            self.ax1[i].set_title(mylab(lab[i]))#,size=8 + 0.7*pscale)
            if i>1:self.ax1[i].set_xlabel(mylab('xx'))#,size=8 +0.7*pscale)#,labelpad=lp)
        
        plt.tight_layout()
        plt.show()
        
        return 

    def build_coeffs(self,sp,ats=None):
        '''
        eta_i=eta^A_i - eta^S_i and idem for rho (rho_i=rho^A_i - rho^S_i)
        From fortran vars in hazel_py.f90:
        !eta_i(1:4)=eta(1:4) - stim(1:4)  
        !rho_i(1:3)=eta(5:7) - stim(5:7)  
        From here: etas=sp.eta[aa,s,:]-sp.stim(aa,s,:) con s =0,1,2,3
                   rhos=idem con s =4,5,6

        Hazel build the total opt coeffs internally for the different 
        synthesis methods but it does not store them in runtime. We do it now, at the end.
        '''
        #from string to hazel.spectrum.Spectrum 
        if type(sp) is not hazel.spectrum.Spectrum:sp=self.spectrum[sp]

        sp.etas=sp.eta[:,0:4,:]-sp.stim[:,0:4,:]
        sp.rhos= np.zeros_like(sp.etas)
        sp.rhos[:,1:4,:]=sp.eta[:,4:7,:]-sp.stim[:,4:7,:]

        codic1={'epsi':sp.eps[:,0,:],'epsq':sp.eps[:,1,:],'epsu':sp.eps[:,2,:],'epsv':sp.stim[:,3,:],
            'etai':sp.etas[:,0,:],'etaq':sp.etas[:,1,:],'etau':sp.etas[:,2,:],'etav':sp.etas[:,3,:],
            'rhoq':sp.rhos[:,0,:],'rhou':sp.rhos[:,1,:],'rhov':sp.rhos[:,2,:]}

        codic2={'eps':sp.eps,'etas':sp.etas,'rhos':sp.rhos}

        return codic1,codic2 

    def plot_coeffs(self,sp,coefs=None,par=None,ats=None,scale=2,figsize=None,tf=4):
        if type(sp) is not hazel.spectrum.Spectrum:sp=self.spectrum[sp]
        
        pscale=self.plotscale
        if scale!=pscale:pscale=scale

        #----------------------------------
        #EDGAR:consider only the atmospheres in sp.
        #self.atms_in_spectrum[sp.name] --->. [['c0'], ['c1','c2']]
        if ats is None:ats=self.atms_in_spectrum[sp.name]#ats=self.atmospheres

        labs=[]
        #get name of atmospheres in sp
        for n, order in enumerate(self.atms_in_spectrum[sp.name] ): #n run layers along the ray
            for k, atm_name in enumerate(order):  #k runs subpixels of topologies c1+c2                                                  
                #at=self.atmospheres[atm]
                labs.append(atm_name)
        lines=[]
        #----------------------------------
        
        cd,cd2=self.build_coeffs(sp) #set sp.etas and sp.rhos
        #cds={**cd, **cd2} #merge the two dictionaries
        
        #font = {'family' : 'normal','weight' : 'bold', 'size'   : 22}
        #matplotlib.rc('font', **font)
        plt.rcParams.update({'font.size': pscale*tf})

        if coefs is None:#default    
            if figsize is not None:fs=figsize
            else:fs=(pscale*4,pscale*3)
            
            lab=['epsi','epsq','epsu','epsv','etai','etaq','etau','etav','','rhoq','rhou','rhov']

            alp=[1.,1.,1.] #make plots of MO terms transparent when not used in the calculation 
            if self.apmosekcl[1]==0:alp[2]=0.3

            f, ax = plt.subplots(nrows=3, ncols=4,figsize=fs,label=self.labelf2)
            for cc,coef in enumerate(['eps','etas','rhos']):
                for k,at in enumerate(ats):
                    for sto in range(4):
                        lx, =ax[cc,sto].plot(sp.wavelength_axis,cd2[coef][k,sto,:],alpha=alp[cc]) 
                        ax[cc,sto].set_title(mylab(lab[4*cc+sto]))
                        ax[cc,sto].set_xlabel(mylab('xx'))
        
                        if (sto==0) and (cc ==2):lines.append(lx)
            f.legend(tuple(lines), tuple(labs), loc=(0.1,0.1), bbox_to_anchor=(0.1, 0.3))
        else:
            f, ax = plt.subplots(nrows=1, ncols=len(coefs),figsize=(pscale*len(coefs),pscale),label=self.labelf2)
            for cc,coef in enumerate(coefs):
                alp=1.#make plots of MO terms transparent when not used in the calculation 
                if (self.apmosekcl[1]=='0' and coef[0:3]=='rho'):alp=0.3
                if coef in cd:
                    for k,at in enumerate(ats):
                        ax[cc].plot(sp.wavelength_axis,cd[coef][k,:],alpha=alp)
                        ax[cc].set_title(mylab(coefs[cc]))
                        ax[cc].set_xlabel(mylab('xx'))

                else:
                    raise Exception("Names can only be epsi,epsq,epsu,epsv,etai,etaq,etau,etav,rhoq,rhou,or rhov.")

        plt.tight_layout()
        plt.show()

        return f,ax


    def plot_funcatmos(self,dlims,hz,atmat=None,axs=None,**pkws): 
        '''
        Compares fractional and not fractional polarization.
        This routine is valid when fractional polarization is not implemented in synthesize_spectrum
        **pkws: remaining keyword arguments for plot_PolyFx
        kwargs={var':'mono','method':2}
        'order' here does not play a role to change the atmosphere, only the reference polynomials
        on screen.
        '''
        scale,tf = 4, 2
        pscale=self.plotscale
        if scale!=pscale:pscale=scale
        plt.rcParams.update({'font.size': pscale*tf})

        axs=self.ax5

        rct=[3,3,9]
        if pkws['plotit']!=9:rct=[4,3,10]

        if axs is None:
            f, self.axs = plt.subplots(nrows=rct[0], ncols=rct[1],figsize=(pscale*2,pscale*2.5))  
            axs=self.axs.ravel()
        else:
            for ax in axs:ax.cla()#delete curves in axes in future calls

        #labs=list(dlims.keys())#[r'$P/I_c$',r'$P/I$']
        #selected labs in set_funcatmos :
        #selected=['B1','B2','B3','tau','v','deltav','beta','a','j10','j20f']# ATMAT ORDER
        #idem but with 'beta' at the end 
        labs=['B1','B2','B3','tau','v','deltav','a','j10','j20f','beta'] #PLOT ORDER
        #But after calling self.set_funcatmos from main, dlims is modified to bunch of labels:
        allabs=labs+['ff', 'nbar']

        for i,ax in enumerate(axs):
            if allabs[i]=='deltav':
                hzi,yi=self.fun_minT([hz[0],hz[-1]],dlims[allabs[i]])
                ax.plot(hzi, yi, '-')
            else:    
                var=pkws['var']
                if allabs[i]=='tau':var='exp'
                self.plot_PolyFx(ax,[hz[0],hz[-1]],dlims[allabs[i]],nps=pkws['nps'],var=var,method=pkws['method'])
                
            ax.set_title(mylab(allabs[i]))
            if i >8:ax.set_xlabel(mylab('hz'))

        #just plot the selected parameters that do change and were generated from a polynomial
        #As beta was stored in atmat after deltav, we must move it to plot it at the end as the labels
        if atmat is not None:
            j=[0,1,2,3,4,5,7,8,9,6]#beta at the end and only 10 values because atmat contains 10
            for i in range(rct[2]):#9 or 10
                axs[i].plot(hz, atmat[j[i],:], 'bo',ms=3)

        #plt.gcf().legend(tuple([l1,l2]), tuple(labs))#loc=(0.1,0.1), bbox_to_anchor=(0.1, 0.3))
        plt.tight_layout()
        plt.show()
        
        return axs


    def get_Tpars(self,lam0,vth=None,temp=None,dlamd=None,atw=4.,vmic=0): 
        '''
        Return all parameters related to Doppler broadening: deltav(i.e. vthermal),T, and dnuD or dlamD
        lam0 in angstroms. Velocities in km/s but we tranform them to cm/s for using cgs constants below
        Microturbulent velocity vmic is optional to transform to/from Temperature
        The equations are:
        vth=np.sqrt(2.0*kb*tt/(mp*atw)+vmic*vmic)  
        T=(vth*vth-vmic*vmic)*mp*atw/(2.0*kb)
        dlamd = (lam0/c)*vth 
        nu0=cc/(lam0);        dnud = (nu0/c)*vth       

        Typical call:
        tem,dlamd,vth=get_Tpars(lam0,atw=self.atwsdic['atom'], vth=vth_array )
        '''
        kb=1.380662E-16   ;mp=1.67265E-24  ;cc=3E10   #CGS
        #kb,cc,mp =1.3807E-23, 3E8, 1.6726E-27 #(J/K=kgm2s-2, m/s, kg)       
        if vth is not None:
            #vth,vmic =  vth*1E5,vmic*1E5  #in cm/s
            temp=1E10*(vth*vth-vmic*vmic)*mp*atw/(2.0*kb) #in Kelvin
            dlamd=1E5*vth*(lam0/cc)  #in Angstroms
        else:
            if temp is not None:
                vth=np.sqrt(2.0*kb*temp/(mp*atw)+vmic*vmic*1E10)# in cm/s  
                dlamd=vth*lam0/cc  #in Angstroms (if lam0 in Angstroms)
                vth=vth*1e-5 #output in km/s
            else: #dlamd must enter in same units as lam0(Angstroms)
                vth=1E-5*dlamd*cc/lam0 #in km/s.  
                temp=1E10*(vth*vth-vmic*vmic)*mp*atw/(2.0*kb) #in Kelvin

        #nu0=cc/(lam0*1E-10)
        #dnud=(nu0/cc)*np.sqrt(2*kb*temp/(matom*mp)+ vmic*vmic*1E10)  

        return vth,temp,dlamd

    def fun_minT(self,hzl,dlims,f2=None,d1=2.,z1=500.,hz=None):
        '''Create non-monotonic function in vth mimicking a chromospehric minimum of T
        at (z1,d1)=(500,2) with exact lower boundary value, minimum value, and an upper boundary value 
        determined by d2=d1*f2. Typically f2 >1.0 for a chromospheric rise 
        and the minimum value at minimum point is vth=2.0, i.e. d1=2.

        Typical call:
        hz=np.linspace(hzlims[0],hzlims[-1],Ncells)
        yfun=fun_minT(hz,4.,f2=1.5)

        '''
        #set data. z1 is hardcoded to 500 km
        if hz is None:hz = np.linspace(hzl[0], hzl[-1], num=30)
        z2=hz[-1]

        d0,d2=dlims[0],dlims[-1]
        if f2 is None:f2=d2/d0
        d2=d0*f2*f2 #f2 squared does the trick to fit d2 approx

        #define function
        BminT= lambda x,a,b,c,gam: (a+b*x+c*np.exp(-gam*x))

        gam=0.02#first guess 0.00125-0.003-0.01  
        #constraints at d1
        bb=(d2-d0*np.exp(-gam*z2))/(z2-(1./gam+z1)*(1.-np.exp(-gam*z2)))
        #consraints at d2
        aa=-bb*(1./gam+z1)
        #constraint at d0
        cc= d0-aa
        #correct gamma assuring constraint at d2
        ff=20.#15-100
        gam1=-(1./z2)*(np.log(d2)-np.log(ff*(d0-aa)))

        argum=-(aa+bb*z1)/(d0-aa)
        if argum>0:gam2=-(1./z1)*np.log(argum)
        gam=np.min([gam2,gam1])

        #recalculate
        bb=(d2-d0*np.exp(-gam*z2))/(z2-(1./gam+z1)*(1.-np.exp(-gam*z2)))
        aa=-bb*(1./gam+z1)
        cc= d0-aa

        #calculate seed function with right shape and relative proportions
        yy=BminT(hz,aa,bb,cc,gam)
        #fit it exactly
        yy=(yy-np.min(yy))
        yyprime=yy*(d0-d1)/yy[0] + d1
        

        return hz,yyprime

    def check_method(self,method):
        if (self.verbose >= 1):self.logger.info('Synthesis method : {0}'.format(method))
        if method not in self.methods_list:raise Exception("I do not know that synthesis method. Stopping.")


    def mutates(self,spec,parsdic=None,\
        B1=None,B2=None,B3=None, tau=None,v=None,deltav=None, beta=None,a=None,ff=None,\
        j10=None,j20f=None,nbar=None, \
        apmosekc=None,dcol=None,compare=True,frac=None,fractional=False,ax=None):
        '''
        EDGAR: New routine to create mutations in the synthesis Model, such that one can repeat an experiment just changing
        a few parameters from a previous synthesis. 
        The versatility of reading pars both from parsdic and from keywords, and the checking of 
        the pars as done originally when setting up the model, make this routine cumbersome.

        apmosekc,dcol,ref_frame, hz  topology,los,boundary,, 
        
        Due to the Python behavior, newmo=self is just a reference assignment where both variable names 
        point to the same object, it would only create newmo as a pointer to the object self, 
        without truly performing an indepedent copy. For doing a copy of arrays one has np.copy(). 
        For objects/dictionaries we have deepcopy. 
     
        Parsdic (and any other mutable type,lists or dictionaries) defined as keyword parameters
        will no reset their values between function calls, so having memory of previous  mutates()
        calls from main. If we do not want to use a DTO class or avoid parsdic, then we need to define its
        default to None and make the following check at the beggining of this function. 
        '''
        if parsdic is None:parsdic={}

        if frac is True:fractional=True #abreviated keyword for fractional

        #EDGAR: Substitute next block by next line and check result
        #keywords = dict(locals())#retrieve keywords as dictionary
        #print(keywords)
        #inefficient way of passing pars through keywords: by copying them to parsdic
        #NOTE: apmosekc is introduced always as keyword, but the specific paraemters to which ap-mo-se-kc
        #makes reference can also be specifically introduced via parsdic dictionary.
        #that is why the four long keywords contained in extrapars_list are not considered in the next
        #block to include them in parsdic, because they are already there and never as individual keywords
        if apmosekc is not None:parsdic['apmosekc']=apmosekc
        if dcol is not None:parsdic['dcol']=dcol
        if B1 is not None:parsdic['B1']=B1
        if B2 is not None:parsdic['B2']=B2
        if B3 is not None:parsdic['B3']=B3
        if tau is not None:parsdic['tau']=tau
        if v is not None:parsdic['v']=v
        if deltav is not None:parsdic['deltav']=deltav
        if beta is not None:parsdic['beta']=beta
        if a is not None:parsdic['a']=a
        if ff is not None:parsdic['ff']=ff
        if j10 is not None:parsdic['j10']=j10
        if j20f is not None:parsdic['j20f']=j20f
        if nbar is not None:parsdic['nbar']=nbar


        #----these vars could be extracted out as global to avoid repetition--------
        message="Unknown mutable parameter. The options are:  \n Atompol, MO effects,Stim. emission, Kill coherences, B1, B2, B3 ...{0}"

        '''Identify Model object parameters to be changed. 
        First, general parameters of the Model: apmosekc, dcol and their product apmosekcl'''
        extrapars_list=['Atompol','MO effects','Stim. emission','Kill coherences']
        #parameters of chromosphere objects
        atmpars=['B1','B2','B3','tau','v','deltav','beta','a','ff','j10','j20f','nbar']
        #parameters that can be mutated for now 
        checkdictio={'apmosekc':None,'dcol':None,'Atompol':None,'MO effects':None,'Stim. emission':None,'Kill coherences':None,
            'B1': None, 'B2': None, 'B3': None,'tau':None,'v':None,'deltav':None,'beta':None,
            'a':None,'ff':None,'j10':None,'j20f':None,'nbar':None}
        #------------------------------------------------------------------------------------------------
        extrapars={}
        mutating_keys=[]

        newmo= copy.deepcopy(self) #here we are already copying the spectrum objects inside spectrum
        
        '''kill ghost figure f1 appearing when replicating the model object.
        Closes eventual open figures/axes to avoid the deep copy of the object plot axes.
        otherwise, an ugly copy of the figure in self.ax0 will pop up when calling again a plt.show''' 
        if newmo.f1 is not None:
            newmo.f1.set_label('deletef1')
            newmo.remove_fig('deletef1')#newmo.f1)
        if newmo.f3 is not None:
            newmo.f3.set_label('deletef3')
            newmo.remove_fig('deletef3')#newmo.f1)


        if (len(parsdic)!=0): #if parsdic is not empty (default is {}) we shall ignore mutation keywords!
            #write a list of all possible parameters and check whether the inputs are compatible/valid
            for elem in parsdic:
                if (elem not in checkdictio):raise Exception(message.format(elem))

            apmosekc=self.apmosekc
            dcol=self.dcol
            #build the extrapars dictionary (if its keywords are in parsdic) to call get_apmosekcl()
            for elem in extrapars_list:
                if (elem in parsdic):extrapars[elem]=parsdic[elem]
            #if (len(extrapars)!=0) or ('apmosekc' in parsdic)
            if ('apmosekc' in parsdic):apmosekc=parsdic['apmosekc']
            if ('dcol' in parsdic):dcol=parsdic['dcol']
            if (len(extrapars)!=0) or ('apmosekc' in parsdic) or ('dcol' in parsdic):#update self.apmosekcl
                newmo.apmosekcl=newmo.get_apmosekcl(apmosekc,dcol,extrapars)
                #newmo.apmosekcl=self.get_apmosekcl(apmosekc,dcol,extrapars)

            for key in atmpars:#key=[atm_number,value]#key value. list mutating keys for atmospheres
                if (key in parsdic):mutating_keys.append(key) #ONLY FOR ATM PARS!

        else:#then we shall read mutation pars only from keywords
            if (apmosekc is not None) or (dcol is not None): #if one of them must be updated
                if (apmosekc is None):apmosekc=self.apmosekc #let the other as in original model
                if (dcol is None):dcol=self.dcol              #let the other as in original model
                newmo.apmosekcl=newmo.get_apmosekcl(apmosekc,dcol,{}) #and update + check
                #newmo.apmosekcl=self.get_apmosekcl(apmosekc,dcol,{}) #and update + check

            #if (B1 is not None):mutating_keys.append(B1)

        #--------------------------------------------------------------------------------
        #Close and reopen Hazel seems to be necessary in order to properly detect the
        #updated values of j10,j20f,and nbar 
        newmo.exit_hazel()   
        newmo.ntrans=hazel_code._init(self.atomfile,0) #We initialize pyhazel (and setup self.ntrans) with verbose=0 

        if type(spec) is not hazel.spectrum.Spectrum:spec=self.spectrum[spec]
        '''        
        original spectrum was already duplicated when duplicated de Model object. 
        Here we just get a pointer to it maintaining the name of spectrum 
        '''
        newspec=newmo.spectrum[spec.name] #just for shortening sintaxis
        '''
        Here we decide to leave same name, but we make here explicit how to proceed otherwise.
        If we wish to change the name of the new spectrum in the new object
        we have to change it everywhere (in self.topologies and in atms_in_spectrum):
        '''
        newspecname=spec.name #same name
        #substitute the key spec.name by newspecname in all dictionaries:
        newmo.spectrum[newspecname]= newmo.spectrum.pop(spec.name) 
        newmo.atms_in_spectrum[newspecname]= newmo.atms_in_spectrum.pop(spec.name)
        newmo.topologies[newspecname]= newmo.topologies.pop(spec.name) #{'sp1':'ch1->ch2'}
       
        '''Now we reset spectrum without calling model.add_spectrum again (inefficient).
        We only need to reset the optical coeffs and stokes, and change few pars in 
        spectrum object: LOS, BOUNDARY.'''        

        #spectrum.add_spectrum (NOT model.add_spectrum) to reset spectrum
        wvl=newspec.wavelength_axis
        wvl_lr=newspec.wavelength_axis_lr 
        wvl_range = [np.min(wvl), np.max(wvl)]#used below
        newspec.add_spectrum(newmo.nch, wvl, wvl_lr)#reset stokes, eps, eta, stim, etas, rhos
        '''
        We could directly modify hazelpars in atmopsheres(with this line in
        chromosphere.py: self.atompol,self.magopt,self.stimem,self.nocoh,self.dcol = hazelpars) 
        without updating apmosekcl in model, but we take advantage of update in Model to make input checks
        '''
        if (self.verbose >= 1):self.logger.info('Mutating atmospheric pars...')
        #run over the existing atmospheres of this spectrum topology and set pars
        for n, order in enumerate(newmo.atms_in_spectrum[newspecname] ): #n run layers along the ray
            for k, atm in enumerate(order):  #k runs subpixels of topologies c1+c2                                                  
                at=newmo.atmospheres[atm]

                """
                Activate this spectrum with add_active_line for all existing atmospheres.
                In normal setup, activate_lines is called after adding all atmospheres in topology.
                """         
                at.add_active_line(spectrum=newspec, wvl_range=np.array(wvl_range))

                #SET HAZELPARS: self.atompol,self.magopt,self.stimem,self.nocoh,self.dcol = hazelpars                
                at.atompol,at.magopt,at.stimem,at.nocoh,at.dcol = newmo.apmosekcl 
                #print(at.atompol,at.magopt,at.stimem,at.nocoh,at.dcol)

                #key=[atm_number,value]#key value. Produce mutation updating atm.dna[key]            
                for key in mutating_keys:#print(key,parsdic[key][0],parsdic[key][1])
                    #if in selected layer mutate one layer at a time for every parameter,
                    #but layer can be different among parameters
                    if (parsdic[key][0]==n+k):at.dna[key]=parsdic[key][1] #mutates dna.  #print(key,n+k)
                pars,kwds=at.get_dna() #get updated dna pars of this single atmosphere                     
                at.set_pars(pars,**kwds)#ff=at.atm_dna[8],j10=at.atm_dna[9],j20f=at.atm_dna[10],nbar=at.atm_dna[11]) 
                        
        #--------------------------------------------------------------------------------
        #Synthesize the new spectrum FROM the new model object:
        newmo.synthesize(method=self.methods_dicT[newmo.synmethod],muAllen=newmo.muAllen,obj=self)
        if (self.verbose >= 1):self.logger.info('Spectrum {0} has mutated.'.format(spec.name))
        
        if (compare is True):self.compare_mutation(spec,newspec,fractional=fractional) 
        
        return newmo,newspec


    def compare_mutation(self,sp,newsp,scale=3,tf=2.4,ls='dashed',fractional=False,
        lab=['iic','qic','uic','vic']):

        #all plots in compare_mutation will always be done with the fractional keyword set 
        #in the very first call to mutates and compare_mutation in main to avoid ill comparisons.
        if self.lock_fractional is None:self.lock_fractional=fractional
        else:fractional=self.lock_fractional

        if fractional:lab=['iic','qi','ui','vi']
        else:lab=['iic','qic','uic','vic']            

        #font = {'family' : 'normal','weight' : 'bold', 'size'   : 22}
        #matplotlib.rc('font', **font)
        pscale=self.plotscale
        if scale!=pscale:pscale=scale
        plt.rcParams.update({'font.size': pscale*tf})
        
        if self.ax3 is None:#create and arrange axes, plot old spectra
            self.f3, ax = plt.subplots(2,2,figsize=(pscale*2,pscale*2),label=self.labelf3)#'Mutations (Model.ax2)')  
            self.ax3 = ax.flatten()    

            #to make it alwayss visible , the original first one is plot last one with dashed black line 
            self.ax3[0].plot(sp.wavelength_axis,sp.stokes[0,:],linestyle=ls,color='k')
        
            for i in range(1,4):
                if fractional:self.ax3[i].plot(sp.wavelength_axis,sp.stokes[i,:]/sp.stokes[0,:],linestyle=ls,color='k')
                else:self.ax3[i].plot(sp.wavelength_axis,sp.stokes[i,:],linestyle=ls,color='k')
            
            for i in range(4):
                self.ax3[i].set_title(mylab(lab[i]))#,size=8 + 0.7*pscale)
                if i>1:self.ax3[i].set_xlabel(mylab('xx'))#,size=8 +0.7*pscale)#,labelpad=lp)

        #plot mutated spectra
        self.ax3[0].plot(newsp.wavelength_axis,newsp.stokes[0,:])
        for i in range(1,4):
            if fractional:self.ax3[i].plot(newsp.wavelength_axis,newsp.stokes[i,:]/newsp.stokes[0,:])
            else:self.ax3[i].plot(newsp.wavelength_axis,newsp.stokes[i,:])

        plt.tight_layout()
        plt.show()

    def use_configuration(self, config_dict):
        """
        Use a configuration file

        Parameters
        ----------
        config_dict : dict
            Dictionary containing all the options from the configuration file previously read

        Returns
        -------
        None
        """        
        # Output file
        self.output_file = config_dict['working mode']['output file']

        # Backtracking mode
        if ('backtracking' in config_dict['working mode']):
            self.backtracking = config_dict['working mode']['backtracking']
        else:
            self.backtracking = 'brent'
        if (self.verbose >= 1):
            self.logger.info('Backtracking mode : {0}'.format(self.backtracking))
        
        
        # Working mode  #EDGAR: when using conf file the keyword is still called working mode
        # self.working_mode = config_dict['working mode']['action']  

        # Deal with the spectral regions        
        #tmp = config_dict['spectral regions']
        # EDGAR: Add spectral regions. 
        #We now move this block after adding atmopsheres below
        #for key, value in config_dict['spectral regions'].items():
        #    self.add_spectral(value)  


        # Deal with the atmospheres----------------------------------------------------------------
        tmp = config_dict['atmospheres']
        self.atmospheres = {}

        if (self.verbose >= 1):
            self.logger.info('Adding atmospheres')

        for key, value in tmp.items():
            
            if ('photosphere' in key):
                if (self.verbose >=1):
                    self.logger.info('  - New available photosphere : {0}'.format(value['name']))

                self.add_photosphere(value)
                                                            
            if ('chromosphere' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available chromosphere : {0}'.format(value['name']))

                self.add_chromosphere(value)
                                            
            if ('parametric' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available parametric : {0}'.format(value['name']))

                self.add_parametric(value)

            if ('straylight' in key):
                if (self.verbose >= 1):
                    self.logger.info('  - New available straylight : {0}'.format(value['name']))

                self.add_straylight(value)
        # ----------------------------------------------------------------
       
        # Add spectral regions. EDGAR: it create self.spectrum and add value['topologies'] to self.topologies
        for key, value in config_dict['spectral regions'].items():
            #EDGAR: if Andres likes the new add_spectrum shift to the second choice
            self.add_spectral(value) #option 1: use old routine add_spectral reading config file as before
            #self.add_spectrum(value['name'],config=value) #option 2: use new routine add_spectrum 
        # ----------------------------------------------------------------
        # Set number of cycles if present
        if (self.working_mode == 'inversion'):
            if ('number of cycles' in config_dict['working mode']):
                if (config_dict['working mode']['number of cycles'] != 'None'):
                    self.n_cycles = int(config_dict['working mode']['number of cycles'])
                    if (self.verbose >= 1):
                        self.logger.info('Using {0} cycles'.format(self.n_cycles))

        # Use analytical RFs if possible
        if ('analytical rf if possible' in config_dict['working mode']):
            if (config_dict['working mode']['analytical rf if possible'] != 'None'):
                self.use_analytical_RF_if_possible = hazel.util.tobool(config_dict['working mode']['analytical rf if possible'])
            else:
                self.use_analytical_RF_if_possible = False
        else:
            self.use_analytical_RF_if_possible = False
        if (self.verbose >= 1):
            self.logger.info('Using analytical RFs if possible : {0}'.format(self.use_analytical_RF_if_possible))

        # Set number of maximum iterations
        if ('maximum iterations' in config_dict['working mode']):
            if (config_dict['working mode']['number of cycles'] != 'None'):
                self.max_iterations = int(config_dict['working mode']['maximum iterations'])
            else:
                self.max_iterations = 10
        else:
            self.max_iterations = 10
        if (self.verbose >= 1):
            self.logger.info('Using {0} max. iterations'.format(self.max_iterations))

        # Randomization
        if (self.verbose >= 1):
            if (self.n_randomization == 1):
                self.logger.info('Not using randomizations')
            else:
                self.logger.info('Using a maximum of {0} randomizations'.format(self.n_randomization))

        # Set number of maximum iterations
        if ('relative error' in config_dict['working mode']):
            if (config_dict['working mode']['relative error'] != 'None'):
                self.relative_error = float(config_dict['working mode']['relative error'])
                if (self.verbose >= 1):
                    self.logger.info('Stopping when relative error is below {0}'.format(self.relative_error))
            else:
                self.relative_error = 1e-4
        else:
            self.relative_error = 1e-4

        # Save all cycles
        if ('save all cycles' not in config_dict['working mode']):
            self.save_all_cycles = False
        else:
            self.save_all_cycles = hazel.util.tobool(config_dict['working mode']['save all cycles'])

        if (self.verbose >= 1):
            self.logger.info('Saving all cycles : {0}'.format(self.save_all_cycles))
        

        #EDGAR: setup requires values in self.topologies that were added in add_spectral
        #we want create global optical coeffs in add_spectral together with stokes var
        #to create global coeffs we require n_chromospheres because is set in setup.
        self.setup()#EDGAR:n_chromospheres is set here. 
        

    def setup(self):
        """
        Setup the model for synthesis/inversion. This setup includes adding the topologies, removing unused
        atmospheres, reading the number of cycles for the inversion and some sanity checks

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
              
        # Adding topologies
        if (self.verbose >= 1):
            self.logger.info("Adding topologies") 
        #for value in self.topologies:
        #    self.add_topology(value)
        for specname, value in self.topologies.items():
            self.add_topology(value,specname)

        # Remove unused atmospheres defined in the configuration file and not in the topology
        if (self.verbose >= 1):
            self.logger.info("Removing unused atmospheres")
        self.remove_unused_atmosphere()        

        # Calculate indices for atmospheres, n_chromospheres is set here
        index_chromosphere = 1
        index_photosphere = 1
        self.n_photospheres = 0
        self.n_chromospheres = 0
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                v.index = index_photosphere
                index_photosphere += 1
                self.n_photospheres += 1
            if (v.type == 'chromosphere'):
                v.index = index_chromosphere 
                index_chromosphere += 1  
                self.n_chromospheres += 1
        #EDGAR:index_chromosphere and index_photosphere could be avoided using v.index as counting variable.
        #Also, as the variation of v.index happens all inside the above loop, that parameter endsup being equal
        #to n_chromospheres, hence most of this loop looks unnecessary.

        if (self.verbose >= 1):#EDGAR: print number of Hazel chromospheres/slabs
            self.logger.info('N_chromospheres at setup',self.n_chromospheres)

        # Use analytical RFs if only photospheres are defined
        if (self.n_chromospheres == 0 and self.use_analytical_RF_if_possible):
            self.use_analytical_RF = True
            if (self.verbose >= 1):
                self.logger.info('Using analytical RFs : {0}'.format(self.use_analytical_RF))
        else:
            self.use_analytical_RF = False


        # Check that number of pixels is the same for all atmospheric files if in synthesis mode
        if (self.working_mode == 'synthesis'):
            n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
            all_equal = all(x == n_pixels[0] for x in n_pixels)
            if (not all_equal):
                for k, v in self.atmospheres.items():
                    self.logger.info('{0} -> {1}'.format(k, v.n_pixel))
                raise Exception("Files with model atmospheres do not contain the same number of pixels")
            else:
                if (self.verbose >= 1):
                    self.logger.info('Number of pixels to read : {0}'.format(n_pixels[0]))
                self.n_pixels = n_pixels[0]

        if (self.working_mode == 'inversion'):
            n_pixels = [v.n_pixel for k, v in self.spectrum.items()]
            all_equal = all(x == n_pixels[0] for x in n_pixels)
            if (not all_equal):
                for k, v in self.spectrum.items():
                    self.logger.info('{0} -> {1}'.format(k, v.n_pixel))
                raise Exception("Files with spectral regions do not contain the same number of pixels")
            else:
                if (self.verbose >= 1):
                    self.logger.info('Number of pixels to invert : {0}'.format(n_pixels[0]))
                self.n_pixels = n_pixels[0]
        

        # Check that the number of pixels from all observations (in case of inversion) is the same
        # Check also if they are equal to those of the models
        # n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
        # all_equal = all(x == n_pixels[0] for x in n_pixels)

        # Check that the number of cycles is the same for all atmospheres (in case of inversion)
        if (self.working_mode == 'inversion'):
            cycles = []
            for k, v in self.atmospheres.items():
                for k2, v2 in v.cycles.items():       
                    if (v2 is not None):
                        cycles.append(len(v2))

            all_equal = all(x == cycles[0] for x in cycles)
            if (not all_equal):
                raise Exception("Number of cycles in the nodes of active atmospheres is not always the same")
            else:
                if (self.n_cycles is None):
                    self.n_cycles = cycles[0]

        
        # if (self.working_mode == 'inversion'):
        #     cycles = []
        #     for tmp in ['I', 'Q', 'U', 'V']:
        #         if ( cycles.append
        #     for k, v in self.atmospheres.items():
        #         for k2, v2 in v.cycles.items():                    
        #             cycles.append(len(v2))

        #     all_equal = all(x == cycles[0] for x in cycles)
        #     if (not all_equal):
        #         raise Exception("Number of cycles in the nodes of active atmospheres is not always the same")
        #     else:
        #         if (self.n_cycles is None):
        #             self.n_cycles = cycles[0]

        filename = os.path.join(os.path.dirname(__file__),'data/LINEAS')
        ff = open(filename, 'r')
        self.LINES = ff.readlines()
        ff.close()

        self.init_sir()

        for k, v in self.spectrum.items():            
            v.allocate_info_cycles(n_cycles=self.n_cycles)

        for k, v in self.atmospheres.items():
            v.allocate_info_cycles(n_cycles=self.n_cycles)

        # Count total number of free parameters
        if (self.working_mode == 'inversion'):
            self.n_free_parameters = 0

            for k, v in self.atmospheres.items():
                for k2, v2 in v.cycles.items():
                    if (v2 is not None):                   
                        self.n_free_parameters += max(hazel.util.onlyint(v2[0:self.n_cycles+1]))

            if (self.verbose >= 1):
                self.logger.info('Total number of free parameters in all cycles : {0}'.format(self.n_free_parameters))

        #if (self.verbose >= 1):#EDGAR: print number of Hazel chromospheres/slabs
        #    self.logger.info('N_chromospheres',self.n_chromospheres)

        #if self.plotit:
        self.setup_set_figure() #self.fig and self.ax are created here
        
        #return f,ax
        #return self.n_chromospheres

    def open_output(self):
        self.output_handler = Generic_output_file(self.output_file)        
        self.output_handler.open(self)

    def close_output(self):        
        self.output_handler.close()

    def write_output(self, randomization=0):
        if (self.working_mode == 'synthesis'):
            self.flatten_parameters_to_reference(cycle=0)        
        self.output_handler.write(self, pixel=0, randomization=randomization)


    def add_spectrum(self, name, config=None, wavelength=None, topology=None, los=None, 
        i0fraction=1.0,boundary=None, atom=None, synmethod=None,
        linehazel=None, linesSIR=None, atmos_window=None, instrumental_profile=None,
        wavelength_file=None, wavelength_weight_file=None,observations_file=None, mask_file=None,
        weights_stokes_i=None,weights_stokes_q=None,weights_stokes_u=None,weights_stokes_v=None):
        """
        Similar to add_spectral but with keywords and more compact.
        Programmatically add a spectral region
        """
        #init keywords that are mutable types : before they were directly inittialized to [None]*10
        #in header, which prevents resetting between function calls 
        if weights_stokes_i is None:weights_stokes_i=[None]*10
        if weights_stokes_q is None:weights_stokes_q=[None]*10
        if weights_stokes_u is None:weights_stokes_u=[None]*10
        if weights_stokes_v is None:weights_stokes_v=[None]*10

        value = dict(locals()) #retrieve keywords as dictionary
        value['name'] =name #we add the positional argument "name" to the dictionary

        #EDGAR: this if block can be deleted if we decide to use the same names 
        #for the config file keywords as for the keywords of the subroutine.
        #Map from config file keywords to compact names in this subroutine. 
        #Case insensitive because that was tackled when reading file
        #In config file we should only define the keywords that we want different than None
        #config is old "value" dictionary read from config file and containing all arguments of this subroutine
        if (config is not None):
            #for k, val in config.items(): --> you run over dictionary elements faster with this            
            if ('Wavelength' in config):wavelength=config['Wavelength']
            if ('topology' in config):topology=config['topology']
            if ('LOS' in config):los=config['LOS']
            if ('i0fraction' in config):i0fraction=config['i0fraction']
            if ('Boundary' in config):boundary=config['Boundary']
            if ('atom' in config):atom=config['atom']
            if ('synmethod' in config):atom=config['synmethod']
            if ('line' in config):linehazel=config['line']
            if ('lines' in config):linesSIR=config['lines']
            if ('spectral region' in config):atmos_window=config['spectral region']
            if ('instrumental profile' in config):instrumental_profile=config['instrumental profile']
            if ('wavelength file' in config):wavelength_file=config['wavelength file']
            if ('wavelength weight file' in config):wavelength_weight_file=config['wavelength weight file']
            if ('observations file' in config):observations_file=config['observations file']
            if ('mask file' in config):mask_file=config['mask file']
            if ('weights stokes i' in config):weights_stokes_i=config['weights stokes i']
            if ('weights stokes q' in config):weights_stokes_q=config['weights stokes q']
            if ('weights stokes u' in config):weights_stokes_u=config['weights stokes u']
            if ('weights stokes v' in config):weights_stokes_v=config['weights stokes v']


        if (self.verbose >= 1):            
            self.logger.info('Adding spectral region {0}'.format(name))        
    
        # Wavelength file is not present
        #if (value['wavelength file'] is None):
        if (wavelength_file is None): #the "is" keyword is used to evaluate whether something is None or not
            # If wavelength is defined            
            if (wavelength is not None): #if ('wavelength' in value):
                axis = wavelength
                wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))                
                wvl_lr = None
                if (self.verbose >= 1):
                    self.logger.info('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
            else:
                raise Exception('Wavelength range is not defined. Please, use "Wavelength" or "Wavelength file"')
        else:
            # If both observed and synthetic wavelength points are given
            if (wavelength is not None):
                axis = wavelength
                if (len(axis) != 3):
                    raise Exception("Wavelength range is not given in the format: lower, upper, steps")
                wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
                if (self.verbose >= 1):
                    self.logger.info('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
                    self.logger.info('  - Reading wavelength axis from {0}'.format(wavelength_file))
                wvl_lr = np.loadtxt(self.root + wavelength_file)
            else:
                if (self.verbose >= 1):
                    self.logger.info('  - Reading wavelength axis from {0}'.format(wavelength_file))
                wvl = np.loadtxt(self.root + wavelength_file)
                wvl_lr = None
                
        if (wavelength_weight_file is None):
            if (self.verbose >= 1 and self.working_mode == 'inversion'):
                self.logger.info('  - Setting all wavelength weights to 1')
            weights = np.ones((4,len(wvl)))
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Reading wavelength weights from {0}'.format(wavelength_weight_file))
            weights = np.loadtxt(self.root + wavelength_weight_file, skiprows=1).T

        # Observations file not present
        if (observations_file is None):
            if (self.working_mode == 'inversion'):
                raise Exception("Inversion mode without observations is not allowed.")            
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using observations from {0}'.format(observations_file))
            
        if (self.verbose >= 1):
            if (mask_file is None):            
                self.logger.info('  - No mask for pixels')
            else:
                self.logger.info('  - Using mask from {0}'.format(mask_file))

        if (self.verbose >= 1):
            if (instrumental_profile is None): 
                self.logger.info('  - No instrumental profile')
            else:
                self.logger.info('  - Instrumental profile : {0}'.format(instrumental_profile))

        if (los is None):
            if (self.working_mode == 'synthesis'):
                raise Exception("You need to provide the LOS for spectral region {0}".format(name))
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using LOS {0}'.format(los))
            los = np.array(los).astype('float64')
        #---------------------------------------------
        if (self.verbose >= 1):
            self.logger.info('  - Using I0fraction = {0} for normalization in spectral region {1}'.format(i0fraction,name))
        '''
        EDGAR: the value of boundary that enters in the call to spectrum below  
        must be Ibackground(physical units)/I0Allen (all input and output Stokes shall always be
        relative to the Allen I0 continuum). Then, the value given to the boundary keyword from outside 
        must be such that is normalized to the Allen continuum. Otherwise, we have to multiply by 
        i0fraction to assure that quantity.
        In general this i0fraction keyword shall be ignored and assumed always to be 1.0, thus concentrating
        all the definition of the boundary condition in the boundary keyword. 

        But if the boundary keyword value was defined in main program 
        as Ibackground(physical units)/Icont_wing (i.e. with the first wing value of I equal to 1.0), 
        then we have to multiply by i0fraction keyword, which should then be different than 1.0 to describe
        a background continuum that is different from Allen. 
        (By the definition of i0fraction we have that Icont_wing is I0fraction*I0Allen).
        
        Example: for boundary=1, the boundary intensity given to Hazel synthesize routine should be
        Iboundary = 1 * I0fraction*I0Allen(lambda). 
        Then, in this section we first multiply boundary * I0fraction, 
        and before calling hazel we shall add units multiplying by I0Allen.
        '''
        if (boundary is None):
            if (self.verbose >= 1):self.logger.info('  - Using default boundary conditions [1,0,0,0] in spectral region {0} or read from file. Check carefully!'.format(name))
            boundary = i0fraction*np.array([1.0,0.0,0.0,0.0])  
            self.normalization = 'on-disk'
        else:
            if (self.verbose >= 1):
                if (np.ndim(boundary[0])==0):#boundary elements are scalars [1.0,0.0,0.0,0.0]
                    self.logger.info('  - Using constant boundary conditions {0}'.format(boundary))
                    if (boundary[0] == 0.0):self.logger.info('  - Using off-limb normalization (peak intensity)')          
                else:#the user already introduced float 64 arrays with spectral dependences for I.
                    self.logger.info('  - Using spectral profiles in boundary conditions')
                    if (boundary[0,0] == 0.0):self.logger.info('  - Using off-limb normalization (peak intensity)')          
            boundary = i0fraction*np.array(boundary).astype('float64')#gives array([1.0,0.0,0.0,0.0]) or array of (4,Nwavelength) 
        
        #---------------------------------------------
        stokes_weights = []
        for st in [weights_stokes_i,weights_stokes_q,weights_stokes_u,weights_stokes_v]:
            tmp = hazel.util.tofloat(st)
            tmp = [i if i is not None else 1.0 for i in tmp]
            stokes_weights.append(tmp)
        stokes_weights = np.array(stokes_weights)
        
        #EDGAR: atom, line_to_index and line keywords moved to add_spectral
        if (atom is not None) and (atom in self.atomsdic):
            self.atom=atom#self.atom can be deleted because is not used anywhere else
            self.line_to_index=self.atomsdic[atom]
        else:
            raise Exception('Atom is not specified or not in the database. Please, define a valid atom.')
    
        if (self.verbose >= 1):self.logger.info('Atom added.')
        
        '''
        EDGAR:The default initialization of self.synmethod was done in the init above. 
        Below is the general set up of the synthesis method to use for synthesizing this spectrum.
        In principle here we assume one single method per spectrum, but one could later 
        change the method during runtime with the keyword method in synthesize(), 
        such that a single spectrum could be synthesize with different methods in different 
        layers during the transfer. Thats is why we do not associate a synmethod to the spectrum object,
        although if necessary we could do it below generating a list of methods in spectrum.
        We work with string names for the user, but internally hazel work with the numbers and the spectrum
        list is made with numbers to shorten because we can have many layers in Hazel2.
        '''
        if (synmethod is not None) and (synmethod != self.methods_dicT[self.synmethod] ):
            self.check_method(synmethod) #synmethod is a string with name, self.synmethod is the number
            self.synmethod=self.methods_dicS[synmethod] #update self.synmethod with the number


        #EDGAR: line for Hazel chromospheres and for SIR photosphere
        #it seems lines for SIR read in add_photosphere were wrong because they were introduced programatically
        #with the field atm['spectral lines'] in add_photosphere, but there was no such a field defined anywhere 
        lineH, lineS = '', ''
        if (linehazel is not None) and (linehazel in self.atomsdic[atom]):
            lineH=linehazel #e.g. '10830'.  Lines for activating in Hazel atmos
            if (self.verbose >= 1):self.logger.info("    * Adding HAZEL line : {0}".format(lineH))
        else:
            if (linesSIR is not None):#same for SIR lines #SIR photospheres have NOT been checked
                lineS = [int(k) for k in list(linesSIR)] #we moved this from add_photosphere
                if (self.verbose >= 1):self.logger.info("    * Adding SIR line : {0}".format(lineS))
            else:
                raise Exception('Line is not specified or not in the database. Please, define a valid line.')

        #Count chromospheres for defining optical coeffs containers, now that all atmospheres have been added
        self.nch=0  #n_chromospheres=0    
        #careful:is this the number of chromospheres added or the number associated to spectrum??
        for k, atm in self.atmospheres.items():            
            if (atm.type == 'chromosphere'):self.nch += 1 #should be equal to self.n_chromospheres.

        if (self.verbose >= 1):self.logger.info('N_chromospheres before setup',self.nch)


        #initialize here the optical coefficient containers with self.nch dimension:
        self.spectrum[name] = Spectrum(wvl=wvl, weights=weights, observed_file=observations_file, 
            name=name, stokes_weights=stokes_weights, los=los, boundary=boundary, 
            mask_file=mask_file, instrumental_profile=instrumental_profile, 
            root=self.root, wvl_lr=wvl_lr,lti=self.line_to_index,lineHazel=lineH,lineSIR=lineS,
            n_chromo=self.nch, synmethod=self.synmethod)

        #EDGAR: update spectrum object with the multiplets for later accesing it from synthesize at chromosphere.py
        self.spectrum[name].multiplets = self.multipletsdic[atom] 
        #ntrans needed to define length of nbar,omega, and j10.
        self.spectrum[name].ntrans = self.ntrans #len(self.multipletsdic[atom])

        #--EDGAR---------------------------------------------------------------
        #we are here defining the wavelength window for all atmospheres associated to this spectral region
        if (atmos_window is not None):#EDGAR: if not in dictionary, then take the one of current spectral region
            wvl_range = [float(k) for k in atmos_window]
        else:
            wvl_range = [np.min(self.spectrum[name].wavelength_axis), np.max(self.spectrum[name].wavelength_axis)]

        #self.topologies.append(topology)#'ph1->ch1+ch2'
        self.topologies[name]=topology# anade una entrada del tipo {'sp1':'ch1->ch2'}
        
        """
        Activate this spectrum with add_active_line for all existing atmospheres.
        Part of this routine was previously inside every add_atmosphere routine.
        Now all spectral and atmospheric actions and routines are disentangled. 
        Activate_lines is now called after adding all atmospheres in topology.
        """
        if (self.verbose >= 1):self.logger.info('Activating lines in atmospheres',self.nch)
        for k, atm in self.atmospheres.items():            
            atm.add_active_line(spectrum=self.spectrum[name], wvl_range=np.array(wvl_range))
                        

        return self.spectrum[name]


    def check_key(dictio,keyword,default):
        if (keyword not in dictio):
            dictio[keyword] = default
        elif (dictio[keyword] == 'None'):
            dictio[keyword] = default
        return dictio

    def add_spectral(self, spectral):
        """
        Programmatically add a spectral region  

        EDGAR:In future versions this subroutine could be deleted in exchange to add_spectrum, which shorter
        and more efficient. Then, remember to call add_spectrum when reading experiment from file.

        Parameters
        ----------
        name: string name of the dictionary
        spectral : dict
            Dictionary containing the following data
            'Name', 'Wavelength', 'Topology', 'Weights Stokes', 'Wavelength file', 'Wavelength weight file',
            'Observations file', 'Mask file','i0fraction','Boundary', 'Synmethod'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        value = hazel.util.lower_dict_keys(spectral)

        #just add the name of the spectrum to the dictionary to keep a fully self-explained dictionary
        #by working directly with name instead of with value['name'] the following code would be more rea
        #value['name']=name 
    
        if (self.verbose >= 1):            
            self.logger.info('Adding spectral region {0}'.format(value['name']))        

        #----EDGAR:much shorter way of checking default values------------- 
        defaultdic={'wavelength file':None,'wavelength weight file':None,'observations file':None,
        'stokes weights':None,'mask file':None,'los':None, 'boundary':None,'i0fraction':1.0,
        'instrumental profile':None,'synmethod':None}#by default i0fraction must be 1.0

        for key in defaultdic:value=check_key(value,key,defaultdic[key])
        #-----------------------------------------------------------

        for tmp in ['i', 'q', 'u', 'v']:
            if ('weights stokes {0}'.format(tmp) not in value):
                value['weights stokes {0}'.format(tmp)] = [None]*10
            elif (value['weights stokes {0}'.format(tmp)] == 'None'):
                value['weights stokes {0}'.format(tmp)] = [None]*10


        # Wavelength file is not present
        if (value['wavelength file'] is None):

            # If the wavelength is defined            
            if ('wavelength' in value):
                axis = value['wavelength']
                wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))                
                wvl_lr = None
                if (self.verbose >= 1):
                    self.logger.info('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
            else:
                raise Exception('Wavelength range is not defined. Please, use "Wavelength" or "Wavelength file"')
        else:
            # If both observed and synthetic wavelength points are given
            if ('wavelength' in value):
                axis = value['wavelength']
                if (len(axis) != 3):
                    raise Exception("Wavelength range is not given in the format: lower, upper, steps")
                wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
                if (self.verbose >= 1):
                    self.logger.info('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
                    self.logger.info('  - Reading wavelength axis from {0}'.format(value['wavelength file']))
                wvl_lr = np.loadtxt(self.root + value['wavelength file'])
            else:
                if (self.verbose >= 1):
                    self.logger.info('  - Reading wavelength axis from {0}'.format(value['wavelength file']))
                wvl = np.loadtxt(self.root + value['wavelength file'])
                wvl_lr = None
                
        if (value['wavelength weight file'] is None):
            if (self.verbose >= 1 and self.working_mode == 'inversion'):
                self.logger.info('  - Setting all wavelength weights to 1')
            weights = np.ones((4,len(wvl)))
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Reading wavelength weights from {0}'.format(value['wavelength weight file']))
            weights = np.loadtxt(self.root + value['wavelength weight file'], skiprows=1).T

        # Observations file not present
        if (value['observations file'] is None):
            if (self.working_mode == 'inversion'):
                raise Exception("Inversion mode without observations is not allowed.")            
            obs_file = None
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using observations from {0}'.format(value['observations file']))
            obs_file = value['observations file']

        if (value['mask file'] is None):            
            mask_file = None
            if (self.verbose >= 1):
                self.logger.info('  - No mask for pixels')
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Using mask from {0}'.format(value['mask file']))
            mask_file = value['mask file']

        if (value['instrumental profile'] is None):
            if (self.verbose >= 1):
                self.logger.info('  - No instrumental profile')
        else:
            if (self.verbose >= 1):
                self.logger.info('  - Instrumental profile : {0}'.format(value['instrumental profile']))

        # if (value['straylight file'] is None):
        #     if (self.verbose >= 1):
        #         self.logger.info('  - Not using straylight')
        #     stray_file = None
        # else:
        #     if (self.verbose >= 1):
        #         self.logger.info('  - Using straylight from {0}'.format(value['straylight file']))
        #     stray_file = value['straylight file']

        if (value['los'] is None):
            if (self.working_mode == 'synthesis'):
                raise Exception("You need to provide the LOS for spectral region {0}".format(value['name']))
            los = None
        else:
            los = np.array(value['los']).astype('float64')
            if (self.verbose >= 1):
                self.logger.info('  - Using LOS {0}'.format(value['los']))

        #this block is adapted from add_spectrum but not checked because add_spectral is deprecated
        #---------------------- 
        if (self.verbose >= 1):
            self.logger.info('  - Using I0fraction = {0} for normalization in spectral region {1}'.format(value['i0fraction'],value['name']))

        

        if (value['boundary'] is None):
            if (self.verbose >= 1):self.logger.info('  - Using default boundary conditions [1,0,0,0] in spectral region {0} or read from file. Check carefully!'.format(value['name']))
            boundary = value['i0fraction']*np.array([1.0,0.0,0.0,0.0])  
            self.normalization = 'on-disk'
        else:
            if (self.verbose >= 1):
                if (np.ndim(value['boundary'][0])==0):#boundary elements are scalars [1.0,0.0,0.0,0.0]
                    self.logger.info('  - Using constant boundary conditions {0}'.format(value['boundary']))
                    if (value['boundary'][0] == 0.0):self.logger.info('  - Using off-limb normalization (peak intensity)')          
                else:#the user already introduced float 64 arrays with spectral dependences for I.
                    self.logger.info('  - Using spectral profiles in boundary conditions')
                    if (value['boundary'][0,0] == 0.0):self.logger.info('  - Using off-limb normalization (peak intensity)')          
            boundary = value['i0fraction']*np.array(value['boundary']).astype('float64')#gives array([1.0,0.0,0.0,0.0]) or array of (4,Nwavelength) 

        if (value['synmethod'] is not None) and (value['synmethod'] != self.methods_dicT[self.synmethod] ):
            self.check_method(value['synmethod']) #synmethod is a string with name, self.synmethod is the number
            self.synmethod=self.methods_dicS[value['synmethod']] #update self.synmethod with the number

        #----------------------

        stokes_weights = []
        for st in ['i', 'q', 'u', 'v']:
            tmp = hazel.util.tofloat(value['weights stokes {0}'.format(st)])
            tmp = [i if i is not None else 1.0 for i in tmp]
            stokes_weights.append(tmp)
        
        stokes_weights = np.array(stokes_weights)
        

        #EDGAR: atom, line_to_index and line keywords moved to add_spectral
        if ('atom' in value) and (value['atom']in self.atomsdic):
            self.atom=value['atom']#self.atom can be deleted because is not used anywhere else
            self.line_to_index=self.atomsdic[value['atom']]
        else:
            raise Exception('Atom is not specified or not in the database. Please, define a valid atom.')
    
        if (self.verbose >= 1):self.logger.info('Atom added.')
        
        #EDGAR: line for Hazel chromospheres and for SIR photosphere
        #it seems lines for SIR read in add_photosphere were wrong because they were introduced programatically
        #with the field atm['spectral lines'] in add_photosphere, but there was no such a field defined anywhere 
        lineH, lineS = '', ''
        if ('linehazel' in value) and (value['linehazel'] in self.atomsdic[value['atom']]):
            lineH=value['linehazel'] #e.g. '10830'.  Lines for activating in Hazel atmos
            if (self.verbose >= 1):self.logger.info("    * Adding HAZEL line : {0}".format(lineH))
        else:
            if ('linesir' in value):#same for SIR lines #SIR photospheres have NOT been checked
                lineS = [int(k) for k in list(value['linesir'])] #we moved this from add_photosphere
                if (self.verbose >= 1):self.logger.info("    * Adding SIR line : {0}".format(lineS))
            else:
                raise Exception('Line is not specified or not in the database. Please, define a valid line.')

        self.nch=0  #n_chromospheres=0    
        for k, atm in self.atmospheres.items():            
            if (atm.type == 'chromosphere'):self.nch += 1 #should be equal to self.n_chromospheres.
        
        if (self.verbose >= 1):self.logger.info('N_chromospheres before setup',self.nch)


        #EDGAR:inside here there is the add_spectrum routine setting up self.spectrum['spx'].wavelength_axis used below 
        self.spectrum[value['name']] = Spectrum(wvl=wvl, weights=weights, observed_file=obs_file, 
            name=value['name'], stokes_weights=stokes_weights, los=los, boundary=boundary, 
            mask_file=mask_file, instrumental_profile=value['instrumental profile'], 
            root=self.root, wvl_lr=wvl_lr,lti=self.line_to_index,lineHazel=lineH,lineSIR=lineS,
            n_chromo=self.nch, synmethod=self.synmethod)
        #we send line_to_index and lines to be activated to this Spectrum object for accessing them later from everywhere
        #we could then remove self.line_to_index from here
        #and we could now remove all what has to do with spectrum from add_chromosphere and add_photosphere
        #so that add_spectral and add_atmos can be invoked in any order


        #--EDGAR---------------------------------------------------------------
        #Update spectrum object with the multiplets for later accesing it from synthesize at chromosphere.py
        self.spectrum[name].multiplets = self.multipletsdic[atom] 
        #ntrans needed to define length of nbar,omega, and j10.
        self.spectrum[name].ntrans = self.ntrans #len(self.multipletsdic[atom])

        #'atmos window'  is the old 'wavelength' keyword of add_chromosphere
        #we are here defining the wavelength window for all atmospheres associated to this spectral region
        value=check_key(value,'atmos window',None)

        if (value[keyw] is not None):#EDGAR: if not in dictionary, then take the one of current atm['spectral region']
            wvl_range = [float(k) for k in value[keyw]]
        else:
            wvl_range = [np.min(self.spectrum[value['name']].wavelength_axis), np.max(self.spectrum[value['name']].wavelength_axis)]

        topo=value['topology']
        #self.topologies.append(topo)#'ph1->ch1+ch2'  #if topologies is a list
        self.topologies[name]=topo# anade una entrada del tipo {'sp1':'ch1->ch2'}  #now topology is dictionary
        
        #this gives just string names:
        #list_of_lists=[k.split('+') for k in topo.split('->')] #[['ph1'], ['ch1', 'ch2']]
        #list_of_atms=[item for sublist in list_of_lists for item in sublist]#['ph1', 'ch1', 'ch2']
    
        """
        Activate this spectrum with add_active_line for all existing atmospheres.
        Part of this routine was previously inside every add_atmosphere routine.
        Now all spectral and atmospheric actions and routines are disentangled. 
        Activate_lines is now called after adding all atmospheres in topology but before passing pars.
        """
        if (self.verbose >= 1):self.logger.info('Activating lines in atmospheres',self.nch)
        for k, atm in self.atmospheres.items():            
            atm.add_active_line(spectrum=self.spectrum[value['name']], wvl_range=np.array(wvl_range))


        #--------------------------------------------------------------------

    def add_photosphere(self, atmosphere):
        """
        Programmatically add a photosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = SIR_atmosphere(working_mode=self.working_mode, name=atm['name'], verbose=self.verbose)
        

        # If NLTE is available because PyTorch and PyTorch Geom are available
        # check whether the line is needed in NLTE or not
        if self.nlte_available:
            if ('nlte' not in atm):
                self.atmospheres[atm['name']].nlte = False
            else:
                self.atmospheres[atm['name']].nlte = hazel.util.tobool(atm['nlte'])
                if (self.verbose >= 1):
                    self.logger.info("    * Line in NLTE if available")
        else:
            self.atmospheres[atm['name']].nlte = False
        
        """  EDGAR:   now you can DELETE this block   
        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        #lines = [int(k) for k in list(atm['spectral lines'])] #we want lines to be set in add_spectral
        #self.atmospheres[atm['name']].add_active_line(lines=lines, spectrum=self.spectrum[atm['spectral region']], 
        #    wvl_range=np.array(wvl_range), verbose=self.verbose) #EDGAR ,  OLD version
                    
        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range), verbose=self.verbose)
        """

        #EDGAR:more efficient set up of reference frame
        refkey=''
        self.atmospheres[atm['name']].reference_frame = 'line-of-sight'#always for photospheres
        if ('reference frame' in atm):refkey='reference frame'
        if ('ref frame' in atm):refkey='ref frame' #short keyword alias
        if ('vertical' in atm[refkey]):
                raise Exception('Magnetic fields in photospheres are always in the line-of-sight reference frame.')

        #if ('reference frame' in atm):
        #    if ('line-of-sight' in atm['reference frame'] or 'LOS' in atm['reference frame']):
        #        self.atmospheres[atm['name']].reference_frame = 'line-of-sight'
        #    if ('vertical' in atm['reference frame']):
        #        raise Exception('Magnetic fields in photospheres are always in the line-of-sight reference frame.')
        #else:
        #    self.atmospheres[atm['name']].reference_frame = 'line-of-sight'

        if (self.verbose >= 1):
            self.logger.info("    * Magnetic field reference frame : {0}".format(self.atmospheres[atm['name']].reference_frame))

        if (self.atmospheres[atm['name']].graphnet_nlte is not None):
            self.set_nlte(True)

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():            
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v                                            
        
        if ('reference atmospheric model' in atm):
            my_file = Path(self.root + atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(self.root + atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()

        if ('nodes' in atm):            
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)                

        if ('temperature change to recompute departure coefficients' in atm):
            self.atmospheres[atm['name']].t_change_departure = float(atm['temperature change to recompute departure coefficients'])
        else:
            self.atmospheres[atm['name']].t_change_departure = 0.0


    def add_chromosphere(self, atmosphere):
        """
        Programmatically add a chromosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength',
            ...
            magnetic field reference frame
            coordinates for magnetic field parameters
            ...
            'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        self.atmospheres[atm['name']]         #None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Hazel_atmosphere(working_mode=self.working_mode, \
        name=atm['name'],ntrans=self.ntrans,hazelpars=self.apmosekcl)#EDGAR:,atom=atm['atom'])

        #------------------------------------------------------------------------
        """EDGAR now you can DELETE THIS BLOCK
        #the only thing of this block that cannot be known from add_spectrum is atm['wavelength']
        #so we need to store it associated to every atmosphere. This is a bit of nonsense because
        #although we would define a spectral window associted to each atmosphere, it is also true
        #that they would all be the same for a same radiative transfer / spectrum topology.
        #so we only have to define a single spectral window ('wavelength') associated to all
        #atmospheres of a same spectral region (or of a same spectrum as we shall rename it).
        #and this have to be done from add_spectrum, although it would be done for each atmosphere,
        #together with the activation of the line.
        
        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):#EDGAR: if not in dictionary, then take the one of current atm['spectral region']
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        #all this could be done in add_spectral (after adding all atmospheres of topology but before passing 
        #parameters) in a for for every atmopshere
        #self.atmospheres[atm['name']].add_active_line(line=atm['line'], spectrum=self.spectrum[atm['spectral region']], 
        #    wvl_range=np.array(wvl_range)) #OLD
        #now line is taken from spectrum, not from atmosphere, we could instead read it from model object adding key line=self.line
        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))
        """
        #------------------------------------------------------------------------

        #EDGAR: more efficient and flexible reference frame set up
        refkey=''
        self.atmospheres[atm['name']].reference_frame = 'vertical'#default
        if ('reference frame' in atm):refkey='reference frame'
        if ('ref frame' in atm):refkey='ref frame' #short keyword alias
        if (refkey != ''):#desired reference frame has been specified
            if (atm[refkey] == 'line-of-sight' or atm[refkey] == 'LOS'):
                self.atmospheres[atm['name']].reference_frame = 'line-of-sight'
            elif (atm[refkey] != 'vertical'):
                raise Exception('Error: wrong specification of reference frame.')

        if (self.verbose >= 1):
            self.logger.info("    * Adding line : {0}".format(atm['line']))
            self.logger.info("    * Magnetic field reference frame : {0}".format(self.atmospheres[atm['name']].reference_frame))

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v
        #EDGAR : now 'coordinates for magnetic field vector' is 'coordB' 
        #but remember to write it here in lower case (see above)
        if ('coordb' in atm):
            if (atm['coordb'] == 'cartesian'):
                self.atmospheres[atm['name']].coordinates_B = 'cartesian'
            if (atm['coordb'] == 'spherical'):
                self.atmospheres[atm['name']].coordinates_B = 'spherical'
        else:
            self.atmospheres[atm['name']].coordinates_B = 'spherical' #default

        self.atmospheres[atm['name']].select_coordinate_system()

        if (self.verbose >= 1):            
            self.logger.info("    * Magnetic field coordinates system : {0}".format(self.atmospheres[atm['name']].coordinates_B))            


        if ('reference atmospheric model' in atm):
            my_file = Path(self.root + atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(self.root + atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        self.atmospheres[atm['name']].height = float(atm['height'])

        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)
        #EDGAR: return directly the dict entry with the name of the chromo that has been just added
        return self.atmospheres[atm['name']] 

    def add_chrom(self, atmosphere):#EDGAR: alias to add_chromosphere
        return self.add_chromosphere(atmosphere)

    def add_Nchroms(self,tags,ckey,hz=None):
        '''
        Add N chromospheres with a list of tags/names(e.g.['c1','c2']), 
        a dictionary of common paraemters (ckeys), and 
        a list of specific heigths as optional keyword.
        '''
        if (hz is None):hz=[0.0]*len(tags)
        #chout=[]
        for kk,ch in enumerate(tags):
            #chout.append(self.add_chromosphere({'name': ch,'height': hz[kk],**ckey}))
            self.chromospheres.append(self.add_chromosphere({'name': ch,'height': hz[kk],**ckey}))

        #return chout, tags
        return self.chromospheres, tags

    def add_funcatmos(self,Ncells,ckey,hzlims=None,hz=None,topo=''):
        '''
        Creates and add a full chromosphere made of N elemental pieces/slabs/cells
        and making certain parameters to vary according to given P-order polinomials.
        There is a dictionary of common parameters (ckeys), and 
        a list of specific heigths as optional keyword. For the moment, this function
        shall just be a small improvement to add_Nchroms.

        Explanation: Hazel2 allows for a versatile serial or parallel combination 
        of slabs that can be associated to different spectral lines, heights, filling 
        factors and topologies for perform radiative trasnfer on them. On the other side,
        other RT codes work directly to stratified fully discretized atmospheres with many 
        points. The function add_funcatmos here is in between these two approaches, 
        pretending to mimick a full atmosphere with several points from serially concatenating
        Hazel atmospheres but yet mantaining a control on the functional variations
        of its physical parameters. The goal is to simplify its creation process because
        when many chromospheres are added as in add_Nchroms, the qualitative properties
        of each of them (like its labels, or reference frame) become irrelevant or equal
        for the whole piece. The next step would be to directly add a realistic atmosphere
        or convert one to a toy full chromosphere as those here built for Hazel. But for that
        it we shall need to work  directly with temperature, density, etc.

        '''
        #create list of tags/names for each cell (['c1','c2',...])
        #actual value of hz only matters in relation with tau. If tau is left defined
        #for each cell as incoming parameter we are also setting the height scale
        #indirectly with dtau=eta_I*dz. Hence we will need to wait for eta_I before returning 
        #a meaningful value of hz scales
        if (hz is None):hz=[0.0]*Ncells 
        self.hzlims=hzlims

        #we must return topology string used later for add spectrum in the main program
        tags=[]
        if topo == '':
            for kk in range(Ncells):
                tags.append('c'+str(kk)) #['c0','c1',...,'c_{N-1}']
                topo=topo+tags[kk]+'->'
            topo=topo[:-2]
        else:#from topo extract tags
            tags=topo.rsplit(sep='->') #['c0','c1',...,'c_{N-1}']
        
        #chout=[]
        for kk in range(Ncells):             
            #chout.append(self.add_chromosphere({'name': tags[kk],'height': hz[kk],**ckey}))
            self.chromospheres.append(self.add_chromosphere({'name': tags[kk],'height': hz[kk],**ckey}))

        #return chout, topo
        return self.chromospheres, topo


    def check_B_vals(self,B1,B2,B3):
        '''
        Check whether magnetic field parameters set in keywords dictionary
        correspond to coordinates (cartesian or spherical) set in the model coordB keyword
        and whether the values are in their corresponding physical ranges
        '''
        anyk=list(self.atmospheres.keys())[0]
        #if (self.coordinates_B == 'spherical'): --> when routine is placed in chromosphere.py 
        if (self.atmospheres[anyk].coordinates_B == 'cartesian'):
            #in cartesian all magnetic field comps goes between -inf(say -3kG) and inf(3kG)
            if not (0.<= np.abs(B1) <= self.limB) or not (0.<= np.abs(B2)<= self.limB) or not (0.<= np.abs(B3)<= self.limB):
                raise Exception('Possible ERROR: values of magnetic field cartesian components not in range?')
        
        if (self.atmospheres[anyk].coordinates_B == 'spherical'):
            #in spherical B\in[0,limB], thetaB\in[0,180], chiB\in[-360,360]
            #chiB is 0,180 but the user has the freedom to choose negative angles and even 
            #angles larger than 360 but that would be unusual
            if not (0.<= np.abs(B1) <=self.limB) or not (0.<= B2 <=180.0) or not (0.0<=np.abs(B3)<=360.0):
                raise Exception('Possible ERROR: values of magnetic field cartesian components not in range?')

        return True


    def fix_point_polyfit_fx(self,n, x, y, xf, yf) :
        '''Solves a system of equations that allow o determine the parameters of a polynomial 
        that fit points (x,y) approximatelly passing exactly through points (xf,yf).
        At the end return the resulting  polynomial'''
        mat = np.empty((n + 1 + len(xf),) * 2)
        vec = np.empty((n + 1 + len(xf),))
        x_n = x**np.arange(2 * n + 1)[:, None]
        yx_n = np.sum(x_n[:n + 1] * y, axis=1)
        x_n = np.sum(x_n, axis=1)
        idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
        mat[:n + 1, :n + 1] = np.take(x_n, idx)
        xf_n = xf**np.arange(n + 1)[:, None]
        mat[:n + 1, n + 1:] = xf_n / 2
        mat[n + 1:, :n + 1] = xf_n.T
        mat[n + 1:, n + 1:] = 0
        vec[:n + 1] = yx_n
        vec[n + 1:] = yf
        params = np.linalg.solve(mat, vec)
        
        sel_pars=params[:n + 1] 

        return np.polynomial.Polynomial(sel_pars)


    def get_yps(self,xps,ypl,sortit=1,method=2,nps=2):
        '''creates few random nps points in y contained between limits ypl 
        creates random polynomial and few points mapping it''' 
        if (np.abs(sortit) != 1): sortit=1#sortit can only be 1 or -1
        if method==1:
            order=4 #for nps=2, large over bumps if order is not 3 or 4
            #polyfx = np.polynomial.Polynomial(np.random.rand(order + 1)) #simpler option
            polyfx = np.polynomial.Polynomial(np.random.uniform(-2, 8, size=(order + 1,)  ))
            y=polyfx(xps)
            yps=ypl[0]+(ypl[1]-ypl[0])*y/np.max(y)
            if ypl[0]>ypl[-1]:yps=np.sort(yps)[::-1]
        
        if method==2:
            #take some random nps-2 points inside the interval
            ypsr = np.random.uniform(low=ypl[0], high=ypl[1], size=nps-2)
            if ypl[0]>ypl[-1]:ypsr=np.sort(ypsr)[::-1]#sort from smaller to larger and reverse
            yps=[ypl[0]]+list(ypsr)+[ypl[1]] #monotonic series always

        return yps

    def get_exp3points(self,xx,xpl,ypl,nps=3):
        from scipy.optimize import curve_fit
        xps = np.linspace(xpl[0],xpl[1],nps)#few nps points contained between limits xpl 
        yps=self.get_yps(xps,ypl,nps=nps)#,sortit=1,method=2)
        par,cov = curve_fit(exp_3points, xps, yps, p0=np.array([0, -1, 1]), absolute_sigma=True)
        return xps,yps,exp_3points(xx,par[0],par[1],par[2])

    def plot_PolyFx(self,ax,xpl,ypl,nps=2,var='mono',method=2): 
        '''This function just plots some reference polynomials and functions to 
        illustrate possible variations to be assigned to the physical variables
        or to visualize them with respect to the actual variations set.
        We can directly work with this function from main program
        nps=2 is below hardcoded for monotonic method
        '''
        #fig = plt.figure()
        #ax = fig.gca()
        npoints=30
        xx = np.linspace(xpl[0], xpl[-1], num=npoints)

        #array of fixed limiting points, AT LEAST including limiting interval points
        xf, yf = np.array(xpl), np.array(ypl)

        if var == 'non-mono':
            #play with the location of the control points to get different variations
            #creates few nps points contained between limits xpl 
            xps = np.linspace(xpl[0],xpl[1],nps)
            yps=self.get_yps(xps,ypl,method=method,nps=nps)#method 2 is preferred by default

            #creates a polynomial function fitting previous points and passing
            # through given fixed points
            for order in [1,2,3,4]:
                myfx=self.fix_point_polyfit_fx(order, xps , yps, xf, yf)
                ax.plot(xx, myfx(xx), '-')
        if var == 'mono':
            #this method works best with nps=2 to deliver monotonic order-N polyn. functions 
            xps = np.linspace(xpl[0],xpl[1],2)  #here only 2 points to achieve monotonicity
            yps=self.get_yps(xps,ypl,method=2)#method 2 preferred by default
            with warnings.catch_warnings():#avoid printing polyfit  warnings
                warnings.simplefilter("ignore")
                for order in [1,2,3,4]:
                    myfx = np.poly1d(np.polyfit(xps, yps, order))
                    ax.plot(xx, myfx(xx), '-')
        if var == 'mint':#bump mimicking minimum of T
            xps = np.linspace(xpl[0],xpl[1],nps)
            yps=self.get_yps(xps,ypl,method=method,nps=nps)#method 2 is preferred by default
            for order in [1,2,3,4]:
                myfx=self.fix_point_polyfit_fx(order, xps , yps, xf, yf)         
                ax.plot(xx, myfx(xx), '-')

        if var == 'exp':   
            xps,yps=[],[]
            myfx=exp_2points(xx,xpl,ypl)#exponential for tau       
            ax.plot(xx, myfx, '-')
            #xps,yps,myfx=self.get_exp3points(xx,xpl,ypl,nps=3)#ypl is dlims['tau']
            #ax.plot(xx, myfx, '-')

        ax.plot(xps, yps, 'bo')
        ax.plot(xf, yf, 'ro')
        plt.show()

        return #ax


    def PolyFx(self,xx,xpl,ypl,nps=2,order=4,npoints=10,var='mono',method=2):
        '''Get order-N polynomial function connecting points with coords xpl ypl.
        Method=monotonic uses a hardcoded nps=2. Test variations
        with plot_PolyFx function.
        Coeffs given by polyfit are in descending order (x**o to x**0).
        xpl,ypl = [p1[0],p2[0]],[p1[1],p2[1]] --> for two points
        ''' 
        #xx = np.linspace(xpl[0], xpl[-1], num=npoints)

        #array of fixed limiting points, AT LEAST including limiting interval points
        xf, yf = np.array(xpl), np.array(ypl)
        
        if var == 'mono':
            xps = np.linspace(xpl[0],xpl[1],2)#few nps points between limits xpl 
            yps=self.get_yps(xps,ypl,method=2)#method 2 is preferred here

            #needs this to avoid printing polyfit  warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                myfx_eval = np.poly1d(np.polyfit(xps, yps, order))
                myfx=myfx_eval(xx)
        else:
            if var=='tau':myfx=exp_2points(xx,xpl,ypl)#exponential for tau  
            #exp3 not working properly:
            if var=='exp3':xps,yps,myfx=self.get_exp3points(xx,xpl,ypl)#ypl is dlims['tau']
            if var=='deltav':aux,myfx=self.fun_minT(0,ypl,hz=xx)#get minimum of temp function
            if var=='non-mono':#crazy non-monotonic
                xps = np.linspace(xpl[0],xpl[1],nps)#few nps points contained between limits xpl 
                yps=self.get_yps(xps,ypl,method=method,nps=nps)#method 2 is preferred by default
                #creates a polynom func fitting previous points and passing
                # through the given fixed points
                myfx_eval=self.fix_point_polyfit_fx(order, xps , yps, xf, yf)
                myfx=myfx_eval(xx)

        return myfx

    def set_funcatm(self,dlims,hztype='lin',hzlims=None,orders=3,**pkws):
        '''EDGAR:Set the parameters of every atmospheric cell assuming a functional polynomial
        variation between ini and end values given as input parameters in dlims. 
        zhlims was stored in model self.hzlims but here we can overwrite with args '''

        #atm=self.chromospheres 

        if dlims.keys():#dict not empty
            for key in ['ff','nbar','v','beta']:#if pars omited in call,here we set them to default constant values
                if dlims.get(key)==None:dlims[key]=[self.dmm[key][1]]*2  #dlims is now complete always
        else: # dict are empty
            raise Exception("Atmosphere cannot be set. I need a dict of parameters in set_funcatm.")
        
        #check the max and min values in input dict parameters:
        for kk in [0,1]:self.check_B_vals(dlims['B1'][kk],dlims['B2'][kk],dlims['B3'][kk])
        
        for key in ['tau','v','deltav','beta','a','ff','j10','j20f','nbar']:#enumerate only non-magn keys
            if (self.dmm[key][0] <= dlims[key][0] <= self.dmm[key][2]) and (self.dmm[key][0] <= dlims[key][1] <= self.dmm[key][2]):ok=1
            else:raise Exception("Atmosphere cannot be set. Revise ini parameters in set_funcatm.")    
        
        #build functions and values for the parameters 
        #tags in the order expected for building the matrix pars2D in the order expected for set_paramaters
        selected=['B1','B2','B3','tau','v','deltav','beta','a','j10','j20f']
        Ncells=len(self.chromospheres)

        if hzlims is not None:self.hzlims=hzlims #overwrite self.hzlims
        else:hzlims=self.hzlims

        hz = np.linspace(hzlims[0], hzlims[-1], Ncells)

        if hztype=='parab':
            xaux=np.linspace(0,Ncells,Ncells)#yaux is hz : yaux=np.linspace(hzlims[0],hzlims[1],Ncells)
            funX=xaux*xaux
            xnew=np.max(xaux)*funX/np.max(funX)#np.exp(-x)
            hz = np.interp(xnew, xaux, hz)


        if type(orders) is int:orders=[orders]*len(selected) 

        pars2D=np.zeros((len(selected),Ncells))
        
        if (Ncells > 2):    
            for kk,key in enumerate(selected):#create polynomial variations
                if (dlims[key][0]==dlims[key][1]):#constant case
                    pars2D[kk,:]=dlims[key][0]
                    if (orders[kk]>0)&(self.verbose >= 1):warnings.warn("The quantity {0} is being forced to keep constant values.".format(key))
                else:
                    if key=='tau' or key=='deltav':#create exponential or minT functions
                        pars2D[kk,:]=self.PolyFx(hz,hzlims,[dlims[key][0],dlims[key][1]],order=orders[kk],npoints=Ncells,var=key)
                    else:
                        pars2D[kk,:]=self.PolyFx(hz,hzlims,[dlims[key][0],dlims[key][1]],order=orders[kk],npoints=Ncells)
        else:  
            for kk,key in enumerate(selected):
                pars2D[kk,:]=np.array( [dlims[key][0],dlims[key][1] ] )

        for ii in range(Ncells):#set chromospheric cells
            self.chromospheres[ii].set_pars(pars2D[0:8,ii],dlims['ff'][1],j10=pars2D[8,ii],j20f=pars2D[9,ii],nbar=dlims['nbar'][1])
        
        if pkws['plotit']!=0:self.plot_funcatmos(dlims,hz,atmat=pars2D,**pkws)
        #if pkws['plotit']:self.plot_funcatmos(dlims,hz,atmat=pars2D,var=pkws['var'],method=pkws['method'])

        return hz,self.chromospheres,pars2D


    def add_parametric(self, atmosphere):
        """
        Programmatically add a parametric atmosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Wavelength', 'Reference atmospheric model', 'Type',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Parametric_atmosphere(working_mode=self.working_mode)
        
        """EDGAR you can now delete this block
        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))
        """

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

        if ('reference atmospheric model' in atm):
            my_file = Path(self.root + atm['reference atmospheric model'])
            if (not my_file.exists()):
                raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

            self.atmospheres[atm['name']].load_reference_model(self.root + atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v

    
    def add_straylight(self, atmosphere):
        """
        Programmatically add a straylight atmosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Straylight_atmosphere(working_mode=self.working_mode)
        
        """
        if ('wavelength' not in atm):
            atm['wavelength'] = None
        elif (atm['wavelength'] == 'None'):
            atm['wavelength'] = None

        if (atm['wavelength'] is not None):
            wvl_range = [float(k) for k in atm['wavelength']]
        else:
            wvl_range = [np.min(self.spectrum[atm['spectral region']].wavelength_axis), np.max(self.spectrum[atm['spectral region']].wavelength_axis)]

        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))
        """

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)        

        my_file = Path(self.root + atm['reference atmospheric model'])
        if (not my_file.exists()):
            raise FileExistsError("Input file {0} for atmosphere {1} does not exist.".format(my_file, atm['name']))

        if ('reference atmospheric model' in atm):
            self.atmospheres[atm['name']].load_reference_model(self.root + atm['reference atmospheric model'], self.verbose)

            if (self.atmospheres[atm['name']].model_type == '3d'):
                self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        if ('nodes' in atm):
            for k, v in atm['nodes'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):                            
                        self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        for k2, v2 in self.atmospheres[atm['name']].parameters.items():
            self.atmospheres[atm['name']].regularization[k2] = None

        if ('regularization' in atm):
            for k, v in atm['regularization'].items():                
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():                    
                    if (k.lower() == k2.lower()):                        
                        if (v == 'None'):
                            self.atmospheres[atm['name']].regularization[k2] = None
                        else:
                            self.atmospheres[atm['name']].regularization[k2] = v
        

    def remove_unused_atmosphere(self):
        """
        Remove unused atmospheres
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        
        to_remove = []        
        for k, v in self.atmospheres.items():
            if (not v.active):
                to_remove.append(k)
                if (self.verbose >= 1):
                    self.logger.info('  - Atmosphere {0} deleted.'.format(k))
                
        for k in to_remove:
            self.atmospheres.pop(k)
                    
    def init_sir_external(self):
        """
        Initialize SIR for this synthesis
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                f = open('lte.grid', 'w')
                f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
                f.write("           b) The first six characters of the last line                \n")
                f.write("          in the header (if any) must contain the symbol ---       \n")
                f.write("\n")                                                                       
                f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
                f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
                f.write("-----------------------------------------------------------------------\n")

                ind_low = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[0])).argmin()
                ind_top = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[1])).argmin()

                low = v.spectrum.wavelength_axis[ind_low]
                top = v.spectrum.wavelength_axis[ind_top]         # TODO
                delta = (v.spectrum.wavelength_axis[1] - v.spectrum.wavelength_axis[0])

                filename = os.path.join(os.path.dirname(__file__),'data/LINEAS')
                ff = open(filename, 'r')
                flines = ff.readlines()
                ff.close()

                for i in range(len(v.lines)):
                    for l in flines:
                        tmp = l.split()
                        index = int(tmp[0].split('=')[0])
                        if (index == v.lines[0]):
                            wvl = float(tmp[2])                    
                                            
                f.write("{0}            :  {1}, {2}, {3}\n".format(str(v.lines)[1:-1], 1e3*(low-wvl), 1e3*delta, 1e3*(top-wvl)))
                f.close()
                
                v.n_lambda = sir_code.init_externalfile(v.index, filename)

    def init_sir(self):
        """
        Initialize SIR for this synthesis. This version does not make use of any external file, which might be
        not safe when running in MPI mode.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        lines = []
        n_lines = 0

        elements = {'H':1,'HE':2,'LI':3,'BE':4,'B':5,'C':6,'N':7,'O':8,'F':9,'NE':10,
            'NA':11,'MG':12,'AL':13,'SI':14,'P':15,'S':16,'CL':17,'AR':18,'K':19,'CA':20,'SC':21,'TI':22,'V':23,'CR':24,
            'MN':25,'FE':26,'CO':27,'NI':28,'CU':29,'ZN':30,'GA':31,'GE':32,'AS':33,'SE':34,'BR':35,'KR':36,
            'RB':37,'SR':38,'Y':39,'ZR':40,'NB':41,'MO':42,'TC':43,'RU':44,'RH':45,'PD':46,'AG':47,'CD':48,'IN':49,
            'SN':50,'SB':51,'TE':52,'I':53,'XE':54,'CS':55,'BA':56,'LA':57,'CE':58,'PR':59,'ND':60,'PM':61,
            'SM':62,'EU':63,'GD':64,'TB':65,'DY':66,'HO':67,'ER':68,'TM':69,'YB':70,'LU':71,'HF':72,'TA':73,'W':74,
            'RE':75,'OS':76,'IR':77,'PT':78,'AU':79,'HG':80,'TL':81,'PB':82,'BI':83,'PO':84,'AT':85,'RN':86,
            'FR':87,'RA':88,'AC':89,'TH':90,'PA':91,'U':92}
        states = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}

        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):

                n_lines += 1
                
                ind_low = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[0])).argmin()
                ind_top = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[1])).argmin()

                low = v.spectrum.wavelength_axis[ind_low]
                top = v.spectrum.wavelength_axis[ind_top]         # TODO
                delta = (v.spectrum.wavelength_axis[1] - v.spectrum.wavelength_axis[0])
                
                nblend = len(v.lines)

                lines = np.zeros(len(v.lines), dtype=np.intc)
                atom = np.zeros(len(v.lines), dtype=np.intc)
                istage = np.zeros(len(v.lines), dtype=np.intc)
                wvl = np.zeros(len(v.lines))
                zeff = np.zeros(len(v.lines))
                energy = np.zeros(len(v.lines))
                loggf = np.zeros(len(v.lines))
                mult1 = np.zeros(len(v.lines), dtype=np.intc)
                mult2 = np.zeros(len(v.lines), dtype=np.intc)
                design1 = np.zeros(len(v.lines), dtype=np.intc)
                design2 = np.zeros(len(v.lines), dtype=np.intc)
                tam1 = np.zeros(len(v.lines))
                tam2 = np.zeros(len(v.lines))
                alfa = np.zeros(len(v.lines))
                sigma = np.zeros(len(v.lines))
                
                for i in range(len(v.lines)):            
                    lines[i] = v.lines[i]                    
                    for l in self.LINES:
                        tmp = l.split()
                        index = int(tmp[0].split('=')[0])
                        if (index == v.lines[i]):
                                                        
                            atom[i] = elements[tmp[0].split('=')[1]]
                            istage[i] = tmp[1]
                            wvl[i] = float(tmp[2])
                            zeff[i] = float(tmp[3])
                            energy[i] = float(tmp[4])
                            loggf[i] = float(tmp[5])
                            mult1[i] = int(tmp[6][:-1])
                            mult2[i] = int(tmp[8][:-1])
                            design1[i] = states[tmp[6][-1]]
                            design2[i] = states[tmp[8][-1]]
                            tam1[i] = float(tmp[7].split('-')[0])
                            tam2[i] = float(tmp[9].split('-')[0])
                            if (len(tmp) == 12):
                                alfa[i] = float(tmp[-2])
                                sigma[i] = float(tmp[-1])
                            else:
                                alfa[i] = 0.0
                                sigma[i] = 0.0
                
                lambda0 = 1e3*(low-wvl[0])
                lambda1 = 1e3*(top-wvl[0])
                n_steps = ind_top - ind_low + 1

                v.n_lambda = n_steps
                
                sir_code.init(v.index, nblend, lines, atom, istage, wvl, zeff, energy, loggf,
                    mult1, mult2, design1, design2, tam1, tam2, alfa, sigma, lambda0, lambda1, n_steps)

    #EDGAR: BUG the exit function defined in pyx file had an undersacore. Now it works
    def exit_hazel(self):
        for k, v in self.atmospheres.items():            
            if (v.type == 'chromosphere'):
                hazel_code._exit(v.index) 

    def add_topology(self, atmosphere_order,specname):
        """
        Add a new topology.
        EDGAR: A topology is always associated to a spectrum. Hence these two aspects should be
        related in such a way that we know exactly the atmospheres by the name of the spectrum 
        (i.e. by the name of the spectral region). If this is done, then in synthesize_spectral_region
        we dont need to run over all atmospheres, but only through those linked to the spectrum name.
        For carrying out the radiative transfer, we need then a routine get_transfer_path 
        that returns the list of atmosphere objects (order) in the self.order_atmospheres list
        associated to the spectrum name.
        

        Parameters
        ----------
        topology : str
            Topology
        
        Returns
        -------
        None

        """

        # Transform the order to a list of lists
        if (self.verbose >= 1):
            self.logger.info('  - {0}'.format(atmosphere_order))

        vertical_order = atmosphere_order.split('->')        
        order = []
        for k in vertical_order:
            name = k.strip().replace('(','').replace(')','').split('+')
            name = [k.strip() for k in name]
            
            tmp = []
            for n in name:
                if (n in self.atmospheres):#EDGAR:check that atmosphres in topology were add before
                    tmp.append(n)
                    self.atmospheres[n].active = True
                else:
                    raise Exception("Atmosphere {0} has not been add. Revise the slab names.".format(name))

            order.append(tmp)
        
        order_flat = [item for sublist in order for item in sublist]

        # Check that straylight components, if any, are not at the last position
        for atm in order_flat[:-1]:            
            if (self.atmospheres[atm].type == 'straylight'):
                raise Exception("Straylight components can only be at the last position of a topology.")
        
        self.order_atmospheres.append(order)
        
        self.atms_in_spectrum[specname]=order  #new for making synthesize_spectral_region easier

        # Check that there are no two photospheres linked with ->
        # because they do not make any sense

        n_photospheres_linked = []
        for atmospheres in self.order_atmospheres:
            for order in atmospheres:
                for k, atm in enumerate(order):
                    if (self.atmospheres[atm].type == 'photosphere'):
                        n_photospheres_linked.append(k)
        
        if (len(n_photospheres_linked) != len(set(n_photospheres_linked))):
            raise Exception("There are several photospheres linked with ->. This is not allowed.")
                        
        

    def normalize_ff(self):
        """
        Normalize all filling factors so that they add to one to avoid later problems.
        We use a softmax function to make sure they all add to one and can be unconstrained

        ff_i = exp(x_i) / sum(exp(x_i))

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """

        for atmospheres in self.order_atmospheres:
            for order in atmospheres:

                total_ff = 0.0
                for atm in order:            
                    if (self.atmospheres[atm].type != 'straylight'):
                        if (self.working_mode == 'inversion'):                            
                            self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1]
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                        else:
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)
                        total_ff += ff

                for atm in order:                    
                    if (self.atmospheres[atm].type != 'straylight'):
                        if (self.working_mode == 'inversion'):
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                            self.atmospheres[atm].parameters['ff'] = ff / total_ff
                            self.atmospheres[atm].parameters['ff'] = physical_to_transformed(self.atmospheres[atm].parameters['ff'], self.atmospheres[atm].ranges['ff'][0], self.atmospheres[atm].ranges['ff'][1])
                        else:
                            ff = transformed_to_physical(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)
                            self.atmospheres[atm].parameters['ff'] = ff / total_ff
                            self.atmospheres[atm].parameters['ff'] = physical_to_transformed(self.atmospheres[atm].parameters['ff'], -0.00001, 1.00001)

    def check_filling_factors(self, spectral_region):
        for n, order in enumerate(self.atms_in_spectrum[spectral_region] ): #n run layers along the ray
            if (len(order) > 1):#k runs subpixels of topologies c1+c2                                                  
                count=0
                for k, atm in enumerate(order):count+=self.atmospheres[atm].parameters['ff']
                if (count!=1.0):
                    print("WARNING: Filling factors of layer {0} do not add up to one. Assuming iso-contribution.".format(n))            
                    for k, atm in enumerate(order):self.atmospheres[atm].parameters['ff']=1.0/len(order)

    def synthesize_spectrum(self, spectral_region, method,perturbation=False, stokes=None,stokes_out = None,fractional=False):
        """
        Synthesize all atmospheres of a spectral region and normalize to continuum of quiet Sun at disk center
        Photospheres must be first lower and unique, stray atms should be always last higher
        Stokes and stokes_out are local variables initialized in header (not intended to be inputs!).
        atms_in_spectrum makes unnecessary to check the asp spectral region inside the double loop below
        -----------Parameters:----------
        spectral_region : str.    Spectral region to synthesize
        method: synthesis method for solving the RTE
        perturbation : bool
            True if you are synthesizing with a perturbation. Then, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        fractional: to calculate emergent Stokes profiles as normalized to continuum or as divided by I(lambda)
        --------------------------------
        for a topology c0->c1+c2 the following inside the loops  gives:
        print(self.atms_in_spectrum[spectral_region])  #0 , [['c0'], ['c1','c2']] #all the chain of atmospehres
        print(n,order)        #0,['c0'] and 1,['c1','c2'] #the current layer
        print(k,atm)            #0,c0 / 0,c1 1,c2 (con chX los nombres(strings) de las atms)
        """        
        for n, order in enumerate(self.atms_in_spectrum[spectral_region] ): #n run layers along the ray
            for k, atm in enumerate(order):  #k runs subpixels of topologies c1+c2                                                  
                self.atmospheres[atm].line_to_index=self.line_to_index#update line_to_index in atm/hazel synthesize with that in add_spectral. 
                asp=self.atmospheres[atm].spectrum #for local compact notation

                xbot, xtop = self.atmospheres[atm].wvl_range
                if (n > 0 and k == 0):stokes_out = stokes[:,xbot:xtop] #Update boundary cond. for layers above bottom one      
            
                if (self.atmospheres[atm].type == 'straylight'):#call to synthesize in stray.py 
                    stokes, error = self.atmospheres[atm].synthesize(nlte=self.use_nlte) 
                    if (error == 1):raise 
                    stokes += (1.0 - self.atmospheres[atm].parameters['ff']) * stokes_out                            
                else: #chromospheres and photospheres
                    if (k == 0):#1st subelement of c1+c2 layer, or just the element in single layers 
                        if (self.use_analytical_RF):stokes, self.rf_analytical, error = \
                            self.atmospheres[atm].synthesize(stokes_out, returnRF=True, nlte=self.use_nlte)#always calls a photosphere with RF  
                        else:stokes, asp.eps[n+k,:,:],asp.eta[n+k,:,:],asp.stim[n+k,:,:],error = \
                            self.atmospheres[atm].synthazel(method,stokes=stokes_out, nlte=self.use_nlte)#For single chromospheres
                    else:#if c1+c2, adds up stokes of all sub-pixels  
                        tmp,asp.eps[n+k,:,:],asp.eta[n+k,:,:],asp.stim[n+k,:,:], error = \
                        self.atmospheres[atm].synthazel(method,stokes=stokes_out,nlte=self.use_nlte)
                        stokes += tmp#DOUBT:the sum of contribs of all supixels should be done at the end of the transfer
                    #-------------------------------------------------------------------
        i0=hazel.util.i0_allen(np.mean(asp.wavelength_axis[xbot:xtop]), self.muAllen)  #at mean wavelength
        #i0=hazel.util.i0_allen(asp.wavelength_axis[xbot:xtop], self.muAllen)[None,:] #at each wavelength
        if (self.use_analytical_RF):#EDGAR CAUTION: verify this is ok out of the loop 
            for k, v in self.rf_analytical.items():
                if (k != 'ff'):v /= i0  #EDGAR CAUTION: fractional keyword does not affect here

        if fractional:i0=stokes[0,:] #when fractional, P(lambda)/I(lambda) will be stored in spectrum object
        
        if (perturbation):asp.stokes_perturbed[:,xbot:xtop] = stokes/ i0
        else:asp.stokes[:,xbot:xtop] = stokes/ i0


    def synthesize_spectral_region(self, spectral_region, perturbation=False): #OLD ROUTINE
        """
        EDGAR:This routine is not called anymore: it has been substituted by synthesize_spectrum.

        Synthesize all atmospheres for a single spectral region and normalize to the continuum of the quiet Sun at disk center

        Parameters
        ----------
        spectral_region : str
            Spectral region to synthesize
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """
        stokes = None
        stokes_out = None

        # Loop over all atmospheres        
        for i, atmospheres in enumerate(self.order_atmospheres):
            for n, order in enumerate(atmospheres):
                for k, atm in enumerate(order):                    
                    #EDGAR: updating line_to_index seen by atmosphere and by hazel synthesize
                    #with the self.line_to_index known from add_spectral: 
                    self.atmospheres[atm].line_to_index=self.line_to_index 
                    asp=self.atmospheres[atm].spectrum

                    if (self.atmospheres[atm].spectrum.name == spectral_region):                        
                        ind_low, ind_top = self.atmospheres[atm].wvl_range
                        # Update the boundary condition only for the first atmosphere if several are sharing ff      
                        if (n > 0 and k == 0):#the multiplication with i0 seems to compensate for the division at the end
                            if (perturbation):
                                stokes_out = asp.stokes_perturbed[:, ind_low:ind_top] * hazel.util.i0_allen(asp.wavelength_axis[ind_low:ind_top], 1.0)[None,:]
                            else:
                                stokes_out = asp.stokes[:, ind_low:ind_top] * hazel.util.i0_allen(asp.wavelength_axis[ind_low:ind_top], 1.0)[None,:]
                        
                        #EDGAR: photospheres must be first lower and unique, stray atms should be always last higher,
                        #and I guess chromospheres can be alone and as many as desired
                        if (self.atmospheres[atm].type == 'straylight'):
                            #EDGAR:these are stray synthesis routines(no stokes_out, and are type 'straylight').
                            stokes, error = self.atmospheres[atm].synthesize(nlte=self.use_nlte) #call to synthesize routine in stray.py 
                            if (error == 1):
                                raise 
                            stokes += (1.0 - self.atmospheres[atm].parameters['ff']) * stokes_out                            
                        else: #for types chromosphere and photosphere
                            if (k == 0): #k=0 is the lower atmosphere, either a PHOTOSPHERE or a CHROMOSPHERE in my tests                               
                                if (self.use_analytical_RF):
                                    stokes, self.rf_analytical, error = self.atmospheres[atm].synthesize(stokes_out, returnRF=True, nlte=self.use_nlte)
                                else:
                                    stokes, asp.eps[n+k,:,:],asp.eta[n+k,:,:],asp.stim[n+k,:,:],error = self.atmospheres[atm].synthesize(stokes=stokes_out, nlte=self.use_nlte)
                            else: 
                                tmp,asp.eps[n+k,:,:],asp.eta[n+k,:,:],asp.stim[n+k,:,:], error = self.atmospheres[atm].synthesize(stokes=stokes_out,nlte=self.use_nlte)
                                stokes += tmp  #EDGAR:we are adding current result to previous last one
                                
                        #-------------------------------------------------------------------
                        # Divide by i0: why are we dividing by i0 in every step of the radiative trasnfer??
                        mean_wvl = np.mean(asp.wavelength_axis[ind_low:ind_top])
                        i0 = hazel.util.i0_allen(mean_wvl, 1.0) #print(i0)

                        if (self.use_analytical_RF):
                            for k, v in self.rf_analytical.items():
                                if (k != 'ff'):v /= i0 

                        if (perturbation):asp.stokes_perturbed[:, ind_low:ind_top] = stokes / i0#[None,:]
                        else:asp.stokes[:, ind_low:ind_top] = stokes / i0#[None,:]
    
    def set_nlte(self, option):
        """
        Set calculation of Ca II 8542 A to NLTE

        Parameters
        ----------
        option : bool
            Set to True to use NLTE, False to use LTE
        """
        self.use_nlte = option
        if (self.verbose >= 1):
            self.logger.info('Setting NLTE for Ca II 8542 A to {0}'.format(self.use_nlte))

    def synthesize(self, perturbation=False, method=None,muAllen=1.0,frac=None,fractional=False,obj=None,plot=None,ax=None):
        """
        Synthesize all atmospheres

        Parameters
        ----------
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """
        if frac is True:fractional=frac #abreviated keyword to fractional

        self.muAllen=muAllen #mu where Allen continuum shall be taken for normalizing Stokes output 

        if (method is not None) and (method != self.methods_dicT[self.synmethod]):
            #print(method,self.methods_dicT[self.synmethod])
            print('Changing synthesis method to {0}.'.format(method))
            self.check_method(method)
            self.synmethod=self.methods_dicS[method]#pass from string label to number label and update self

        #EDGAR: WARNING,I think normalize_ff was not working for the synthesis. 
        if (self.working_mode == 'inversion'):
            self.normalize_ff()
            fractional=False #always work with Stokes/Icont in inversions.

        for k, v in self.spectrum.items():#k is name of the spectrum or spectral region
            #EDGAR: checking correct filling factors in composed layers
            #TBD: this kind of check should be done during setup, not in calculations time 
            self.check_filling_factors(k)
            
            #self.synthesize_spectral_region(k, perturbation=perturbation)  #OLD          
            self.synthesize_spectrum(k, self.synmethod,perturbation=perturbation)   #EDGAR NEW
            #we never call synthesize with fractional=True to avoid storing the fractional
            #polarization in spectrum and thus avoid possible mistakes
            #the fractional polarization shall only be shown in plotting

            if (v.normalization == 'off-limb'):
                if (perturbation):
                    v.stokes_perturbed /= np.max(v.stokes_perturbed[0,:])
                else:
                    v.stokes /= np.max(v.stokes[0,:])

            if (v.psf_spectral is not None):                
                for i in range(4):
                    if (perturbation):
                        v.stokes_perturbed[i,:] = scipy.signal.convolve(v.stokes_perturbed[i,:], v.psf_spectral, mode='same', method='auto')
                    else:
                        v.stokes[i,:] = scipy.signal.convolve(v.stokes[i,:], v.psf_spectral, mode='same', method='auto')
            
            if (v.interpolate_to_lr):
                for i in range(4):
                    if (perturbation):                        
                        v.stokes_perturbed_lr[i,:] = np.interp(v.wavelength_axis_lr, v.wavelength_axis, v.stokes_perturbed[i,:])
                    else:                        
                        v.stokes_lr[i,:] = np.interp(v.wavelength_axis_lr, v.wavelength_axis, v.stokes[i,:])                    

            
            if (plot is not None):#plot called inside loops
                if (k == plot):self.plot_stokes(plot,fractional=fractional)
            else:#plot is None because synthesize routine was called without intention of plotting or from mutation
                if obj is None:self.remove_fig(self.labelf1)#self.f1)
                
            
            
        #return ax

    def find_active_parameters(self, cycle):
        """
        Find all active parameters in all active atmospheres in the current cycle

        Parameters
        ----------
        cycle : int
            Cycle to consider
        
        Returns
        -------
        None

        """
        pars = []
        coupled = []
        self.nodes = []
        left = 0
        right = 0
        for atmospheres in self.order_atmospheres:
            for n, order in enumerate(atmospheres):
                for k, atm in enumerate(order):                    
                    for l, par in self.atmospheres[atm].cycles.items():
                        if (par is not None):
                            if (hazel.util.isint(par[cycle])):
                                if (par[cycle] > 0):
                                    
                                    # [Atmosphere name, n_nodes, nodes, value, range]
                                    self.atmospheres[atm].nodes[l] = np.zeros(par[cycle])

                                    self.atmospheres[atm].n_nodes[l] = par[cycle]

                                    right += par[cycle]
                                    
                                    n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                    tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 
                                        'ranges': self.atmospheres[atm].ranges[l], 'delta': self.atmospheres[atm].epsilon[l],
                                        'left': left, 'right': right, 'regularization': self.atmospheres[atm].regularization[l],
                                        'coupled': False}

                                    self.nodes.append(self.atmospheres[atm].nodes[l])
                                        
                                    left = copy.copy(right)
                                                                                                            
                                    pars.append(tmp)

                                else:
                                    self.atmospheres[atm].nodes[l] = 0.0
                                    self.atmospheres[atm].n_nodes[l] = 0
                            else:
                                
                                n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 'coupled': True}

                                coupled.append(tmp)

        self.active_meta = pars
        self.coupled_meta = coupled

        if (not self.nodes):
            raise Exception("No parameters to invert in cycle {0}. Please add them or reduce the number of cycles. ".format(cycle))            

        self.nodes = np.concatenate(self.nodes).ravel()
        
        
    def synthesize_and_compute_rf(self, compute_rf=False, include_jacobian=False):
        """
        Compute response functions for all free parameters according to all active_parameters

        Parameters
        ----------
        compute_rf : bool (optional, default False)
            If True, then compute the response functions. If not, just compute the synthesis.
        
        Returns
        -------
        None

        """
        if (not compute_rf):            
            return

        n_active_pars = len(self.active_meta)

        loop = 0
        loop2 = 0

        self.hessian_regularization = np.zeros(self.n_free_parameters_cycle)
        self.grad_regularization = np.zeros(self.n_free_parameters_cycle)

        # self.use_analytical_RF = False

        for par in self.active_meta:
            nodes = self.nodes[par['left']:par['right']]

            lower = par['ranges'][0]
            upper = par['ranges'][1]

            if (self.verbose >= 4):
                self.logger.info(" * RF to {0} - {1} - nodes={2}".format(par['parameter'], par['atm'], par['n_nodes']))

            if (self.use_analytical_RF):
                for i in range(par['n_nodes']):
                    rf = {}
                    for k, v in self.spectrum.items():

                        # The minus sign comes from the fact that we compute the RF numerically as
                        # (stokes-stokes_perturbed)/delta
                        # rf[k] = -self.atmospheres[par['atm']].rf_analytical[par['parameter']][:,:,i] * jacobian
                        rf[k] = -self.rf_analytical[par['parameter']][:,:,i]
                        rf[k] = rf[k][None, :, :]

                    if (loop == 0):
                        self.response = rf
                    else:
                        for k, v in self.spectrum.items():
                            self.response[k] = np.vstack([self.response[k], rf[k]])

                    if (par['regularization'] is not None):
                        if (par['regularization'][0] == 'l2-value'):
                            alpha = float(par['regularization'][1])
                            lower = par['ranges'][0]
                            upper = par['ranges'][1]
                            value = physical_to_transformed(float(par['regularization'][2]), lower, upper)                            
                            self.grad_regularization[par['left']:par['right']] = 2.0 * alpha * (self.atmospheres[par['atm']].nodes[par['parameter']] - value)
                            self.hessian_regularization[par['left']:par['right']] = 2.0 * alpha

                    loop += 1
                
            else:           
                for i in range(par['n_nodes']):
                    perturbation = np.zeros(par['n_nodes'])
                    if (nodes[i] == 0):
                        perturbation[i] = self.epsilon * par['delta']
                    else:
                        perturbation[i] = self.epsilon * nodes[i]

                    # Perturb this parameter
                    self.atmospheres[par['atm']].nodes[par['parameter']] = nodes + perturbation

                    # Also perturb those parameters that are coupled
                    for par2 in self.coupled_meta:
                        if (par2['coupled'] is True):                            
                            if (par['atm'] == par2['n_nodes'] and par['parameter'] == par2['parameter']):
                                if (self.verbose >= 4):
                                    self.logger.info("   * Coupling RF to {0} - {1}".format(par2['parameter'], par2['atm']))
                                self.atmospheres[par2['atm']].nodes[par2['parameter']] = nodes + perturbation
                                    
                    # Synthesize
                    self.synthesize(perturbation=True)
                                                    
                    # And come back to the original value of the nodes
                    self.atmospheres[par['atm']].nodes[par['parameter']] = nodes

                    for par2 in self.coupled_meta:
                        if (par2['coupled'] is True):
                            if (par['atm'] == par2['n_nodes'] and par['parameter'] == par2['parameter']):
                                self.atmospheres[par2['atm']].nodes[par2['parameter']] = nodes
                    
                    if (include_jacobian):
                        # jacobian =
                        # self.atmospheres[par['atm']].jacobian[par['parameter']]                    
                        jacobian = jacobian_transformation(nodes[i], lower, upper)                    
                    else:
                        jacobian = 1.0

                    rf = {}                    
                    for k, v in self.spectrum.items():
                        if (v.interpolate_to_lr):
                            rf[k] = jacobian * np.expand_dims((v.stokes_lr - v.stokes_perturbed_lr) / perturbation[i], 0)
                        else:
                            rf[k] = jacobian * np.expand_dims((v.stokes - v.stokes_perturbed) / perturbation[i], 0)

                
                    # rf = np.expand_dims((self.spectrum['spec1'].stokes - self.spectrum['spec1'].stokes_perturbed) / perturbation[i], 0)
                    
                    if (loop == 0):                        
                        self.response = rf
                    else:
                        # self.response = np.vstack([self.response, rf])
                        for k, v in self.spectrum.items():
                            self.response[k] = np.vstack([self.response[k], rf[k]])

                    if (par['regularization'] is not None):
                        if (par['regularization'][0] == 'l2-value'):
                            alpha = float(par['regularization'][1])
                            lower = par['ranges'][0]
                            upper = par['ranges'][1]
                            value = physical_to_transformed(float(par['regularization'][2]), lower, upper)                            
                            self.grad_regularization[par['left']:par['right']] = 2.0 * alpha * (self.atmospheres[par['atm']].nodes[par['parameter']] - float(par['regularization'][2]))
                            self.hessian_regularization[par['left']:par['right']] = 2.0 * alpha

                    loop += 1

        #     for i in range(par['n_nodes']):
        #         rf = {}
        #         for k, v in self.spectrum.items():

        #             # The minus sign comes from the fact that we compute the RF numerically as
        #             # (stokes-stokes_perturbed)/delta
        #             # rf[k] = -self.atmospheres[par['atm']].rf_analytical[par['parameter']][:,:,i] * jacobian
        #             rf[k] = -self.rf_analytical[par['parameter']][:,:,i]
        #             rf[k] = rf[k][None, :, :]

        #         if (loop2 == 0):
        #             self.response2 = rf
        #         else:
        #             for k, v in self.spectrum.items():
        #                 self.response2[k] = np.vstack([self.response2[k], rf[k]])

        #         loop2 += 1

        # import matplotlib.pyplot as pl
        # f, ax = pl.subplots(nrows=3, ncols=2, figsize=(9,9))
        # ax = ax.flatten()
        # for i in range(3):
        #     ax[i].plot(self.response['spec1'][i,0,0:60], label='numerical')
        #     ax[i].plot(self.response2['spec1'][i,0,0:60], label='analytical')
        # ax[i].legend()
        # pl.show()
        # breakpoint()

        # # self.response = copy.deepcopy(self.response2)

        # self.use_analytical_RF = True        

                
    def flatten_parameters_to_reference(self, cycle):
        """
        Flatten all current parameters to the reference atmosphere

        Parameters
        ----------
        cycle : int
            Current cycle
        
        Returns
        -------
        None

        """                
        if (self.working_mode == 'inversion'):
            for k, v in self.atmospheres.items():
                v.set_reference(cycle=cycle)
            
        for k, v in self.spectrum.items():            
            v.stokes_cycle[cycle] = copy.deepcopy(v.stokes)
            if (v.interpolate_to_lr):
                v.stokes_lr_cycle[cycle] = copy.deepcopy(v.stokes_lr)

            if (self.working_mode == 'inversion'):                
                v.chi2_cycle[cycle] = copy.deepcopy(v.chi2)                        
                v.bic_cycle[cycle] = copy.deepcopy(self.n_free_parameters * np.log(v.dof) + v.dof * np.log(v.rss))
                v.aic_cycle[cycle] = copy.deepcopy(2.0 * self.n_free_parameters + v.dof * np.log(v.rss))
        
    def set_new_model(self, nodes):
        """
        Set the nodes of the current model to the values passed on the arguments

        Parameters
        ----------
        nodes : float
            Array with the new set of nodes
        
        Returns
        -------
        None

        """

        n_active_pars = len(self.active_meta)  

        # Modify all active parameters
        for par in self.active_meta:
            left = par['left']
            right = par['right']

            self.atmospheres[par['atm']].nodes[par['parameter']] = nodes[left:right]

        # Modify all coupled parameters accordingly
        for par in self.coupled_meta:
            for par2 in self.active_meta:
                if (par2['atm'] == par['n_nodes'] and par2['parameter'] == par['parameter']):
                    
                    left = par2['left']
                    right = par2['right']
                    
                    self.atmospheres[par['atm']].nodes[par['parameter']] = nodes[left:right]
                    self.atmospheres[par['atm']].parameters[par['parameter']] = copy.copy(self.atmospheres[par2['atm']].parameters[par2['parameter']])
                    

    def modified_svd_inverse(self, H, tol=1e-8):
        """
        Compute the inverse of the Hessian matrix using a modified SVD, by thresholding each subpsace separately

        Parameters
        ----------
        H : float
            Hessian matrix

        tol : float
            Tolerance for the singular value of each subspace
        
        Returns
        -------
        None

        """

        try:
            U, w, VT = np.linalg.svd(H, full_matrices=False)
        except np.linalg.LinAlgError:
            U, w, VT = scipy.linalg.svd(H, full_matrices=False, lapack_driver='gesvd')   # This calculation should be more robust but slower

        w_new = np.zeros_like(w)
        
        for par in self.active_meta:            
            left = par['left']
            right = par['right']

            Ui = np.zeros_like(U)
            Ui[:,left:right] = U[:,left:right]

            Gamma_i = np.diagonal(np.diag(w) @ Ui.T @ U).copy()
            
            wmax = np.max(np.abs(Gamma_i))            
            Gamma_i[np.abs(Gamma_i) < tol*wmax] = 0.0

            w_new += Gamma_i

        w_new_inv = np.zeros_like(w)
        ind = np.where(w_new != 0)[0]
        w_new_inv[ind] = 1.0 / w_new[ind]
        
        return U, w_new_inv, VT
        

    def compute_chi2(self, only_chi2=False, weights=None):
        """
        Compute chi2 for all spectral regions

        Parameters
        ----------
        obs : float
            Vector of observations
        only_chi2 : bool
            Control whether the gradient and Hessian is returned
        
        Returns
        -------
        None

        """
        chi2 = 0.0
        rss = 0.0
        n = len(self.nodes)
        dchi2 = np.zeros(n)
        ddchi2 = np.zeros((n,n))

        for k, v in self.spectrum.items():
            if (v.interpolate_to_lr):
                residual = (v.stokes_lr - v.obs)
            else:
                residual = (v.stokes - v.obs)
            
            # Do not use weights. This is used for the computation of errors            
            # if (weights is None):
            weights = (v.stokes_weights[:,self.cycle][:,None] * v.wavelength_weights) * v.factor_chi2            

            chi2 += np.sum(weights * residual**2)

            rss += np.sum(residual**2)
            
            if (not only_chi2):
                response = self.response[k]
                dchi2 += -2.0 * np.sum(weights[None,:,:] * response * residual[None,:,:] , axis=(1,2)) #/ v.dof
                ddchi2 += 2.0 * np.sum(weights[None,None,:,:] * response[None,:,:,:] * response[:,None,:,:] , axis=(2,3)) #/ v.dof

            v.chi2 = chi2
            v.rss = rss
            
        if (not only_chi2):                
            return chi2, dchi2, ddchi2
        else:                
            return chi2

        # if (not only_chi2):            
        #     return chi2, dchi2, ddchi2
        # else:
        #     return chi2

    def compute_uncertainty(self):
        """
        Compute the uncertainty in the parameters at the minimum with the current Hessian

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
    
        #----------------------------
        # Recalculate chi2 without weights
        # chi2 = 0.0
        # for k, v in self.spectrum.items():
        #     residual = (v.stokes - v.obs)

        # weights = v.dof / residual**2

        # # Calculate Hessian
        # self.synthesize_and_compute_rf(compute_rf=True)
        # chi2, dchi2, ddchi2 = self.compute_chi2(weights=weights)
        # hessian = 0.5 * ddchi2


        #----------------------------
        # Recalculate chi2 without weights        

        # Calculate Hessian
        self.synthesize_and_compute_rf(compute_rf=True, include_jacobian=True)
        chi2, dchi2, ddchi2 = self.compute_chi2()
        hessian = 0.5 * ddchi2
        
        U, w_inv, VT = self.modified_svd_inverse(hessian, tol=self.svd_tolerance)
        cov = VT.T.dot(np.diag(w_inv)).dot(U.T)

        # breakpoint()

        for par in self.active_meta:
            left = par['left']
            right = par['right']

            dof = self.atmospheres[par['atm']].spectrum.dof
            rf = scipy.stats.chi2(dof)
            delta = np.sqrt(rf.isf(1.0 - scipy.special.erf(1.0/np.sqrt(2.0))))
            
            rf = scipy.stats.chi2(right-left)
            delta = np.sqrt(rf.isf(1.0 - scipy.special.erf(1.0/np.sqrt(2.0))))

            cov_diagonal = np.abs(np.diagonal(cov[left:right,left:right]))

            # This gives 1sigma error in the transformed domain
            error = np.sqrt(cov_diagonal) * delta

            # Multiply by the Jacobian of the transformation to compute the error in the physical quantities
            error *= jacobian_transformation(self.nodes[left:right], par['ranges'][0], par['ranges'][1])
                    
            self.atmospheres[par['atm']].error[par['parameter']] = error

    def _fun_backtracking(self, log_lambda, dchi2, ddchi2):
        H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
        H += np.diag(10.0**(log_lambda) * np.diag(H))
        gradF = 0.5 * (dchi2 + self.grad_regularization)

        U, w_inv, VT = self.modified_svd_inverse(H, tol=self.svd_tolerance)

        # xnew = xold - H^-1 * grad F
        delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
        
        # Clip the new solution so that the step is resaonable
        new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
    
        self.set_new_model(new_solution)
                    
        self.synthesize_and_compute_rf()
        
        chi2 = self.compute_chi2(only_chi2=True)

        if (self.verbose >= 4):
            self.logger.info('  - Backtracking - lambda: {0:7.5f} - chi2: {1:7.5f}'.format(10.0**log_lambda, chi2))
                    
        return chi2
            
    
    def backtracking_brent(self, dchi2, ddchi2, maxiter=10, bounds=[-3.0,3.0], tol=1e-2):        
        tmp = scipy.optimize.minimize_scalar(self._fun_backtracking, bounds=bounds, args=(dchi2, ddchi2), 
            method='bounded', options={'xatol': tol, 'maxiter': maxiter})

        return 10.0**tmp['x']
        
    def backtracking_parabolic(self, dchi2, ddchi2, direction='down', maxiter=5, lambda_init=1e-3, current_chi2=1e10):
        """
        Do the backtracking to get an optimal value of lambda in the LM algorithm

        Parameters
        ----------
        dchi2 : float
            Gradient of the chi2
        ddchi2 : float
            Second order derivatives with which the Hessian is computed
        direction : str, optional
            Direction on which do the backtracking ('down'/'up' for decreasing/increasing lambda)
        maxiter : int
            Maximum number of iterations
        lambda_init : float
            Initial value of lambda
        current_chi2 : float
            Current best chi2 to compare with those of the backtracking
        
        Returns
        -------
        lambda_opt : float
            Optimal value of lambda found. Bracketed value if bracketing has been possible or just the best value otherwise
        bracketed : bool
            True if the best value has been bracketed
        best_chi2 : float
            Best value of chi2 found
        """
        
        lambdaLM = lambda_init

        chi2_arr = []
        lambdas = []
        sols = []
        keepon = True
        bracketed = False
        loop = 0
        best_chi2 = current_chi2

        while keepon:

            H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
            H += np.diag(lambdaLM * np.diag(H))
            gradF = 0.5 * (dchi2 + self.grad_regularization)

            U, w_inv, VT = self.modified_svd_inverse(H, tol=self.svd_tolerance)

            # xnew = xold - H^-1 * grad F
            delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
            
            # Clip the new solution so that the step is resaonable
            new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
            sols.append(new_solution)
        
            self.set_new_model(new_solution)
                        
            self.synthesize_and_compute_rf()
            
            chi2_arr.append(self.compute_chi2(only_chi2=True))
            
            lambdas.append(lambdaLM)

            if (self.verbose >= 4):
                if (direction == 'down'):
                    self.logger.info('  - Backtracking: {0:2d} - lambda: {1:7.5f} - chi2: {2:7.5f}'.format(loop, lambdaLM, chi2_arr[-1]))
                else:
                    self.logger.info('  * Backtracking: {0:2d} - lambda: {1:7.5f} - chi2: {2:7.5f}'.format(loop, lambdaLM, chi2_arr[-1]))
            
            # If we improve the chi2
            if (chi2_arr[-1] < best_chi2):
                best_chi2 = chi2_arr[-1]

            ind_min = np.argmin(chi2_arr)
            
            if (loop > 1):
                
                # Have we bracketed the minimum
                if (ind_min != 0 and ind_min != len(chi2_arr)-1):
                    keepon = False
                    bracketed = True

            # If lambda < 1e-3, then stop
            if (lambdaLM < 1e-3 or loop > maxiter):
                keepon = False
                min_found = False
            
            if (direction == 'down'):
                lambdaLM /= np.sqrt(10.0)
            else:
                lambdaLM *= np.sqrt(10.0)
            
            loop += 1

        # Parabolic interpolation of the optimal value of lambda
        if (bracketed):
            coeff = np.polyfit(np.log(lambdas[ind_min-1:ind_min+2]), chi2_arr[ind_min-1:ind_min+2], 2)
            lambda_opt = np.exp(-coeff[1] / (2.0*coeff[0]))
        else:
            lambda_opt = lambdas[ind_min]

        return lambda_opt, bracketed, best_chi2, np.min(chi2_arr)

    def randomize(self):
        """
        Randomize all free parameters to lie uniformly in the interval [-2,2] in the transformed
        domain
        """        
        self.nodes = np.random.uniform(low=-2.0, high=2.0, size=self.nodes.shape)

    def invert(self, randomize=False, randomization_ind=None):
        """
        Invert all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        
        first = True

        # Reset reference model to the one loaded from the file
        for k, v in self.atmospheres.items():
            v.reset_reference()

        # Compute normalization factor for the chi^2    
        for k, v in self.spectrum.items():
            v.factor_chi2 = 1.0 / (v.noise**2 * v.dof)        
        
        lambdaLM = 10.0
        lambda_opt = 10.0
        bestchi2 = 1e10
    
        for self.cycle in range(self.n_cycles):            

            if (self.verbose >= 2):
                self.logger.info('-------------')
                if (randomization_ind):
                    self.logger.info('  Cycle {0} - Randomization {1} '.format(self.cycle, randomization_ind))
                else:
                    self.logger.info('  Cycle {0}  '.format(self.cycle))
                
                for k, v in self.spectrum.items():
                    self.logger.info('  Weights for region {0} : SI={1} - SQ={2} - SU={3} - SV={4}'.format(k, v.stokes_weights[0,self.cycle], v.stokes_weights[1,self.cycle],
                        v.stokes_weights[2,self.cycle], v.stokes_weights[3,self.cycle]))
                                

                self.logger.info('-------------')
                

            # Find all active parameters for this cycle and print them in the output
            self.find_active_parameters(self.cycle)

            tmp = [pars['atm'] for pars in self.active_meta]
            tmp = list(set(tmp))

            self.n_free_parameters_cycle = 0

            
            for k, v in self.atmospheres.items():
                if (k in tmp):
                    if (self.verbose >= 3):
                        self.logger.info('Free parameters for {0}'.format(k))
                    for pars in self.active_meta:
                        if (pars['atm'] == k):
                            if (self.verbose >= 3):                                
                                if (pars['coupled'] is False):
                                    if (pars['n_nodes'] == 1):
                                        if (pars['regularization'] is not None):
                                            self.logger.info('  - {0} with {1} node - Regularization -> type:{2}, weight:{3}, value:{4}'.format(pars['parameter'], 
                                                pars['n_nodes'], pars['regularization'][0], pars['regularization'][1], pars['regularization'][2]))
                                        else:
                                            self.logger.info('  - {0} with {1} node - Not regularized'.format(pars['parameter'], pars['n_nodes']))
                                    else:
                                        if (pars['regularization'] is not None):
                                            self.logger.info('  - {0} with {1} nodes - Regularization -> type:{2}, weight:{3}, value:{4}'.format(pars['parameter'], 
                                                pars['n_nodes'], pars['regularization'][0], pars['regularization'][1], pars['regularization'][2]))
                                        else:
                                            self.logger.info('  - {0} with {1} nodes - Not regularized'.format(pars['parameter'], pars['n_nodes']))
                                else:
                                    self.logger.info('  - {0} coupled to {1} variable'.format(pars['parameter'], pars['n_nodes']))
                            if (pars['coupled'] is False):
                                self.n_free_parameters_cycle += pars['n_nodes']
                        
            # Randomize parameters if necessary
            if (randomize):
                self.randomize()

            keepon = True
            iteration = 0

            # Main Levenberg-Marquardt algorithm
            self.synthesize_and_compute_rf(compute_rf=True)
            chi2, dchi2, ddchi2 = self.compute_chi2()       

            while keepon:                                
                
                # Simple parabolic backtracking
                if (self.backtracking == 'parabolic'):
                    lambda_opt, bracketed, best_chi2, backtracking_bestchi2_down = self.backtracking(dchi2, ddchi2, direction='down', maxiter=5, lambda_init=lambdaLM, current_chi2=chi2)

                    backtracking_bestchi2 = copy.copy(backtracking_bestchi2_down)

                    # If solution is not bracketed, then try on the other sense and use the best of the two
                    if (not bracketed):
                        lambda_opt_up, bracketed, best_chi2_up, backtracking_bestchi2_up = self.backtracking(dchi2, ddchi2, direction='up', maxiter=2, lambda_init=lambdaLM)                    

                        if (best_chi2_up < best_chi2):
                            lambda_opt = lambda_opt_up

                        backtracking_bestchi2 = np.min([backtracking_bestchi2, backtracking_bestchi2_up])

                # Bounded Brent backtracking
                if (self.backtracking == 'brent'):                    
                    lambda_opt = self.backtracking_brent(dchi2, ddchi2, maxiter=10, bounds=[-4.0,1.0], tol=1e-2)
                                                
                # if (self.verbose >= 3):
                    # self.logger.info('  * Optimal lambda: {0}'.format(lambda_opt))
                                            
                # If after backtracking the chi2 is larger than the current one, then increase lambda and go to the iteration
                # print(chi2, backtracking_bestchi2)
                # if (chi2 < backtracking_bestchi2 and iteration > 1):
                #     lambdaLM *= 100.0
                #     # print('breaking')
                #     continue


                # Give the final step
                H = 0.5 * (ddchi2 + np.diag(self.hessian_regularization))
                H += np.diag(lambda_opt * np.diag(H))
                gradF = 0.5 * (dchi2 + self.grad_regularization)
                
                U, w_inv, VT = self.modified_svd_inverse(H, tol=self.svd_tolerance)

                # xnew = xold - H^-1 * grad F
                delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)

                # New solution
                # Clip the new solution so that the step is resaonable
                new_solution = self.nodes + np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)
                self.set_new_model(new_solution)

                # Clip the new solution so that the step is resaonable
                self.nodes += np.clip(delta, -self.step_limiter_inversion, self.step_limiter_inversion)

                self.synthesize_and_compute_rf(compute_rf=True)

                chi2, dchi2, ddchi2 = self.compute_chi2()
                

                rel = 2.0 * (chi2 - bestchi2) / (chi2 + bestchi2)

                if (self.verbose > 2):
                    for k, v in self.atmospheres.items():
                        self.logger.info('')
                        self.logger.info('-----------')
                        self.logger.info('{0}'.format(k))
                        self.logger.info('-----------')
                        if (v.type == 'chromosphere'):
                            v.print_parameters(first=first)
                        if (v.type == 'photosphere'):
                            v.print_parameters(first=first)
                        if (v.type == 'parametric'):
                            v.print_parameters(first=first)
                    first = False

                if (self.verbose >= 2):
                    self.logger.info('==============================================================================')
                    self.logger.info('It: {0} - chi2: {1:10.6f} - lambda_opt: {2:10.6f} - rel: {3:10.6f}'.format(iteration, chi2, lambda_opt, np.abs(rel)))
                    self.logger.info('==============================================================================')

                # Increase the optimal by 100 to find again the optimal value
                lambdaLM = 100.0 * lambda_opt

                bestchi2 = copy.copy(chi2)

                if (np.abs(rel) < self.relative_error or iteration > self.max_iterations):
                    keepon = False

                iteration += 1
                                        
            self.set_new_model(self.nodes)

            # Calculate final chi2
            # self.synthesize_and_compute_rf()
            # chi2 = self.compute_chi2(only_chi2=True)
            

            self.compute_uncertainty()
            # if (self.verbose >= 2):
            #     self.atmospheres['ch1'].print_parameters(first=first, error=True)

            self.flatten_parameters_to_reference(self.cycle)

    def _func_grad(self, x):
        """
        Auxiliary functions to use with optimization methods that use gradients
        """
        self.nodes = x
        self.set_new_model(self.nodes)
        self.synthesize_and_compute_rf(compute_rf=True)
        self.chi2, dchi2, _ = self.compute_chi2()
        return self.chi2, dchi2

    def _func_nograd(self, x):
        """
        Auxiliary functions to use with optimization methods that do not use gradients
        """        
        self.nodes = x
        self.set_new_model(self.nodes)
        self.synthesize_and_compute_rf(compute_rf=False)
        self.chi2 = self.compute_chi2(only_chi2=True)
        return self.chi2

    def _callback_general(self, x):
        if (self.verbose >= 2):
            self.logger.info('chi2: {0}'.format(self.chi2))

    def invert_external(self, algorithm, use_jacobian=False, **kwargs):
        """
        Invert all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """

        for k, v in self.spectrum.items():
            v.factor_chi2 = 1.0 / (v.noise**2 * v.dof)
            
        for self.cycle in range(self.n_cycles):
            if (self.verbose >= 2):
                self.logger.info('-------------')
                self.logger.info('  Cycle {0}  '.format(self.cycle))
                
                for k, v in self.spectrum.items():
                    self.logger.info('  Weights for region {0} : SI={1} - SQ={2} - SU={3} - SV={4}'.format(k, v.stokes_weights[0,self.cycle], v.stokes_weights[1,self.cycle],
                        v.stokes_weights[2,self.cycle], v.stokes_weights[3,self.cycle]))

                self.logger.info('-------------')
                

            self.find_active_parameters(self.cycle)

            tmp = [pars['atm'] for pars in self.active_meta]
            tmp = list(set(tmp))

            self.n_free_parameters_cycle = 0
            
            for k, v in self.atmospheres.items():
                if (k in tmp):
                    if (self.verbose >= 3):
                        self.logger.info('Free parameters for {0}'.format(k))
                    for pars in self.active_meta:
                        if (pars['atm'] == k):
                            if (self.verbose >= 3):
                                if (pars['coupled'] is False):
                                    if (pars['n_nodes'] == 1):
                                        self.logger.info('  - {0} with {1} node'.format(pars['parameter'], pars['n_nodes']))
                                    else:
                                        self.logger.info('  - {0} with {1} nodes'.format(pars['parameter'], pars['n_nodes']))
                                else:
                                    self.logger.info('  - {0} coupled to {1} variable'.format(pars['parameter'], pars['n_nodes']))
                            if (pars['coupled'] is False):
                                self.n_free_parameters_cycle += pars['n_nodes']

            n_pars = len(self.nodes)

            if (use_jacobian):
                tmp = algorithm(self._func_grad, self.nodes, jac=True, callback=self._callback_general, **kwargs)
            else:
                tmp = algorithm(self._func_nograd, self.nodes, callback=self._callback_general, **kwargs)

            self._func_grad(tmp['x'])
            self._callback_general(tmp['x'])

            self.set_new_model(tmp['x'])

            self.flatten_parameters_to_reference(self.cycle)
            
    def read_observation(self):
        for k, v in self.spectrum.items():
            v.read_observation(pixel=self.pixel)
            # v.read_straylight(pixel=self.pixel)
