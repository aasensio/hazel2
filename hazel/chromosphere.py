from collections import OrderedDict
import numpy as np
import os
from hazel.atmosphere import General_atmosphere
from hazel.util import i0_allen
from hazel.codes import hazel_code
from hazel.hsra import hsra_continuum
from hazel.io import Generic_hazel_file
from hazel.exceptions import NumericalErrorHazel
import copy


__all__ = ['Hazel_atmosphere']

class Hazel_atmosphere(General_atmosphere):
    def __init__(self, working_mode, name='', ntrans=0, hazelpars=None):
    
        super().__init__('chromosphere', name=name)

        #check/init mutable keywords to avoid memory preservation among function/class calls
        if hazelpars is None:hazelpars=[1,1,1,0,[0.,0.,0.]]

        self.height = 3.0
        self.working_mode = working_mode
        self.ntr=ntrans #taking ntrans from model object
        
        #EDGAR. j10 and j20f(omega) dynamically dimensioned to the number of transitions
        #For now, they are not in dict of parameters because probably will not work in inversions 
        #because the inversion requires these pars to be scalars , not vectors. 
        self.j10 = np.zeros(self.ntr)  #this is j10
        self.j20f = np.ones(self.ntr)  #these are reduction factors, not j20
        self.nbar = np.ones(self.ntr)  #these are reduction factors, not j20
        #EDGAR: this is commented because we need to evaluate whether is possible
        #to do inversions with radiation field anisotropies for every transition
        #self.parameters['j10'],self.parameters['j20f'] = np.zeros(self.ntr) ,np.ones(self.ntr)  
        self.atompol,self.magopt,self.stimem,self.nocoh,self.dcol = hazelpars  #extract the pars for Hazel 
        
        self.dna= {}#None #init memory of pars for doing mutations faster (it is updated in set_pars below)


        #EDGAR:much more compact way of initializing parameters
        keyparlist=['Bx','By','Bz','B','thB','phiB','tau','v','deltav','beta','a','ff'] #,'j10', 'j20f'
        unitslist=['G', 'G', 'G',  'G','deg','deg',  '', 'km/s','km/s', '',   '',  '']  #,'%'  , '' 

        for k,key in enumerate(keyparlist):
            self.units[key]=unitslist[k]        
            self.nodes[key]=0.0                
            self.nodes_location[key]= None
            self.error[key] = 0.0        
            self.n_nodes[key] = 0                
            self.ranges[key] = None        
            self.cycles[key] = None
            self.epsilon[key] = 0.01                
            self.jacobian[key] = 1.0
            self.regularization[key] = None                
            self.parameters[key]=0.0

        #only parameters initialized to something different than 0.0
        self.parameters['beta'],self.parameters['tau'] = 1.0, 1.0
        self.parameters['deltav'] = 8.0
        self.parameters['ff'] = np.log(1.0)


    def select_coordinate_system(self):
        #EDGAR: if you choose one coordB, this pops the other one out
        if (self.coordinates_B == 'cartesian'):
            labels = ['B', 'thB', 'phiB']
        if (self.coordinates_B == 'spherical'):
            labels = ['Bx', 'By', 'Bz']
        
        for label in labels:            
            self.parameters.pop(label)
            # self.units.pop(label)
            # self.nodes.pop(label)
            # self.error.pop(label)
            # self.n_nodes.pop(label)
            # self.ranges.pop(label)
            # self.cycles.pop(label)
            # self.jacobian.pop(label)
            # self.regularization.pop(label)
    

    def add_active_line(self, spectrum, wvl_range):#EDGAR def add_active_line(self, line, spectrum, wvl_range):

        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------        
        lines : str
            Line to activate: ['10830','5876']
        spectrum : Spectrum
            Spectrum object
        wvl_range : float
            Vector containing wavelength range over which to synthesize this line
        
        Returns
        -------
        None
    
        EDGAR: number of transitions is already in atom%ntran
        """      
        #EDGAR:removed this block from here once line_to_index is defined in add_spectral
        #if (self.atom == 'helium')  :
        #    self.line_to_index = {'10830': 1, '3888': 2, '7065': 3,'5876': 4}
        #if (self.atom == 'sodium')  :
        #    self.line_to_index = {'5895': 1, '5889': 2}
            

        self.active_line =spectrum.lineHazel # line
        self.wavelength_range = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top+1]


    def cartesian_to_spherical(self,Bx,By,Bz):
        # Transform to spherical components in the vertical reference frame which are those used in Hazel
        B = np.sqrt(Bx**2 + By**2 + Bz**2)
        if (B == 0):
            thetaB = 0.0
        else:
            thetaB = 180.0 / np.pi * np.arccos(Bz / B)
        phiB = 180.0 / np.pi * np.arctan2(By, Bx)
    
        return B, thetaB,phiB

    def los_to_vertical(self,Bx_los,By_los,Bz_los):
        Bx_vert = Bx_los * self.spectrum.mu + Bz_los * np.sqrt(1.0 - self.spectrum.mu**2)
        By_vert = By_los
        Bz_vert = -Bx_los * np.sqrt(1.0 - self.spectrum.mu**2) + Bz_los * self.spectrum.mu

        return Bx_vert,By_vert,Bz_vert


    def spherical_to_cartesian(self,B,thetaB,phiB):
        Bx = B * np.sin(thetaB * np.pi / 180.0) * np.cos(phiB * np.pi / 180.0)
        By = B * np.sin(thetaB * np.pi / 180.0) * np.sin(phiB * np.pi / 180.0)
        Bz = B * np.cos(thetaB * np.pi / 180.0)

        return Bx,By,Bz        


    def get_B_Hazel(self,B1,B2,B3):#pars[0],pars[1],pars[2]
        '''
        Transform magnetic field parameters to vertical reference frame 
        in spherical coordinates, which is the Hazel working system.
        Output provides both cartesian and spherical coordinates
        in case both are needed (I would leave just spherical).
        '''
        if (self.coordinates_B == 'cartesian'):
            # Cartesian - LOS : just carry out the rotation in cartesian geometry and transform to spherical
            if (self.reference_frame == 'line-of-sight'):
                Bx,By,Bz=self.los_to_vertical(B1,B2,B3) #here Bi would be cartesian Bx,By,Bz in LOS
            else:# Cartesian - vertical : do nothing
                Bx,By,Bz=B1,B2,B3

            # Transform to spherical components in the vertical reference frame which are those used in Hazel
            B,thB,phiB=self.cartesian_to_spherical(Bx,By,Bz) #from cartesian vertical to spherical vertical

        if (self.coordinates_B == 'spherical'):
            # Spherical - vertical by default: do nothing but get cartesians .
            # we assume vertical but result is overwritten in next block if it is not like that
            Bx,By,Bz=self.spherical_to_cartesian(B1,B2,B3) #delete this line if Bx,By,Bz not needed in the code 
            B,thB,phiB= B1,B2,B3

            # Spherical - LOS : transform to cartesian, do rotation to vertical and come back to spherical
            if (self.reference_frame == 'line-of-sight'):
                Bx_los,By_los,Bz_los=self.spherical_to_cartesian(B,thB,phiB)
                Bx,By,Bz=self.los_to_vertical(Bx_los,By_los,Bz_los)
                # To spherical comps in vertical frame as required in Hazel
                B,thB,phiB=self.cartesian_to_spherical(Bx,By,Bz)

        return B,thB,phiB,Bx,By,Bz #Bx,By,Bz can be deleted if not needed in the code 

    def get_dna(self):
        #extract dna of previous experiment in a way that is readable by set_pars directly
        #j10, j20f, and nbar are stored before always as a list of ntrans elements
                
        #extract dna of atmosphere of previous experiment
        pars=[val[1] for val in self.dna.items()][0:8]
        kwds={'ff':self.dna['ff'],'j10':self.dna['j10'], 'j20f':self.dna['j20f'],'nbar': self.dna['nbar']}

        return pars,kwds


    def set_parameters(self, pars, ff=1.0,j10=None,j20f=None,nbar=None,m=None):
        """
        Set the parameters of this model chromosphere

        Parameters 
        ----------
        pars : list of float
            This list contains the following parameters in order: Bx, By, Bz, tau, v, delta, beta, a

        ff : float optional keyword
            Filling factor

        #EDGAR: 
        j10: array of doubles optional keyword, in percentage units. It is a vector with Ntransitions length
            but in this header we init it to have 10 transitions because it cannot be predefined with 
            self.ntr in the definition of a subroutine.
        Returns
        -------
        None
        """
        self.parameters['B'],self.parameters['thB'],self.parameters['phiB'], \
        self.parameters['Bx'],self.parameters['By'],self.parameters['Bz']= \
        self.get_B_Hazel(pars[0],pars[1],pars[2]) 

        self.parameters['tau'] = pars[3]
        self.parameters['v'] = pars[4]
        self.parameters['deltav'] = pars[5]
        self.parameters['beta'] = pars[6]
        self.parameters['a'] = pars[7]
        self.parameters['ff'] = ff
        
        #EDGAR: j10 and j20f can be introduced as a float that will here be broadcasted to the number of 
        #transitions, or directly as a list with Ntransitions (self.ntr) length.
        if (j10 is not None):#then overwrite its default zero value set above with self.j10=np.zeros(self.ntr)
            if (type(j10) is not list):
                self.j10 = np.zeros(self.ntr)+ j10 #if not a list, assumed a number that sets array to j10 keyword
            else:#is a list
                if len(j10)!=self.ntr:raise Exception('ERROR: j10 should have {0} elements'.format(self.ntr) )
                self.j10 = np.array(j10) #from list to array of doubles              

        if (j20f is not None):#then overwrite its default 1.0 value set above
            if (type(j20f) is not list):
                self.j20f = np.zeros(self.ntr)+ j20f #if not a list, assumed a number that sets array to j20f keyword
            else:#is a list
                if len(j20f)!=self.ntr:raise Exception('ERROR: j20f should have {0} elements'.format(self.ntr) )
                self.j20f = np.array(j20f) #from list to array of doubles              

        if (nbar is not None):#then overwrite its default 1.0 value set above
            if (type(nbar) is not list):
                self.nbar = np.zeros(self.ntr)+ nbar #if not a list, assumed a number that sets array to j10 keyword
            else:#is a list
                if len(nbar)!=self.ntr:raise Exception('ERROR: nbar should have {0} elements'.format(self.ntr) )
                self.nbar = np.array(nbar) #from list to array of doubles              


        #CREATES DNA combining pars with ff,j10j20f,nbar
        #AS LIST: add extra pars as lists to the pars list--> dna=pars.append(ff,list(self.j10),list(self.j20f),list(self.nbar)) 
        #AS DICTIONARY with keys 'B1','B2','B3','tau','v','deltav','beta','a','ff', 'j10','j20f','nbar' 
        #we are storing the last extrapars as introduced by user,  without the above modifications
        self.dna['B1']=pars[0]
        self.dna['B2']=pars[1]
        self.dna['B3']=pars[2]
        self.dna['tau']=pars[3]
        self.dna['v']=pars[4]         #dna memory exposed to / required for future mutations
        self.dna['deltav']=pars[5]
        self.dna['beta']=pars[6]
        self.dna['a']=pars[7]
        self.dna['ff']=ff
        self.dna['j10']=j10
        self.dna['j20f']=j20f
        self.dna['nbar']=nbar


        #EDGAR: I think this is not being used because ranges are defined to none above.
        # Check that parameters are inside borders by clipping inside the interval with a border of 1e-8
        if (self.working_mode == 'inversion'):
            for k, v in self.parameters.items():
                if (self.ranges[k] is not None):
                    self.parameters[k] = np.clip(v, self.ranges[k][0] + 1e-8, self.ranges[k][1] - 1e-8)
        

    #alias to set_parameters.:
    def set_pars(self, pars,ff=1.0,j10=None,j20f=None,nbar=None,m=None):
        return self.set_parameters(pars,ff,j10=j10,j20f=j20f,nbar=nbar,m=m)

    def load_reference_model(self, model_file, verbose):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            verbosity flag

        Returns
        -------
            None
        """
        extension = os.path.splitext(model_file)[1][1:]
        
        if (extension == '1d'):
            if (verbose >= 1):
                self.logger.info('    * Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            self.model_filename = model_file          
        
        if (extension == 'h5'):
            if (verbose >= 1):
                self.logger.info('    * Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            

        self.model_handler = Generic_hazel_file(model_file)
        self.model_handler.open()
        out, ff = self.model_handler.read(pixel=0)
        self.model_handler.close()

        self.set_parameters(out, ff)

        self.init_reference(check_borders=True)

    def nodes_to_model(self):
        """
        Transform from nodes to model
        
        Parameters
        ----------
        None
                                
        Returns
        -------
        None
        """         
        for k, v in self.nodes.items():            
            # if (self.n_nodes[k] > 0):                
            self.parameters[k] = self.reference[k] + np.squeeze(self.nodes[k])
            # else:                
                # self.parameters[k] = self.reference[k]            
                            
    def print_parameters(self, first=False, error=False,pre='8.3f'): #EDGAR: now you can modify the precision at once here
        if (self.coordinates_B == 'cartesian'):
            self.logger.info("     {0}        {1}        {2}        {3}       {4}       {5}      {6}      {7}".format('Bx', 'By', 'Bz', 'tau', 'v', 'deltav', 'beta', 'a'))
            self.logger.info("{0:{pp}}  {1:{pp}}  {2:{pp}}  {3:{pp}}  {4:{pp}}  {5:{pp}}  {6:{pp}}  {7:{pp}}".format(self.parameters['Bx'], self.parameters['By'], self.parameters['Bz'], self.parameters['tau'], 
                self.parameters['v'], self.parameters['deltav'], self.parameters['beta'], self.parameters['a'],pp=pre))
            
            if (error):            
                self.logger.info("{0:{pp}}  {1:{pp}}  {2:{pp}}  {3:{pp}}  {4:{pp}}  {5:{pp}}  {6:{pp}}  {7:{pp}}".format(self.error['Bx'], self.error['By'], self.error['Bz'], self.error['tau'], 
                self.error['v'], self.error['deltav'], self.error['beta'], self.error['a'],pp=pre))
        
        if (self.coordinates_B == 'spherical'):
            self.logger.info("     {0}        {1}        {2}        {3}       {4}       {5}      {6}      {7}".format('B', 'thB', 'phiB', 'tau', 'v', 'deltav', 'beta', 'a'))
            self.logger.info("{0:{pp}}  {1:{pp}}  {2:{pp}}  {3:{pp}}  {4:{pp}}  {5:{pp}}  {6:{pp}}".format(self.parameters['B'], self.parameters['thB'], self.parameters['phiB'], self.parameters['tau'], 
                self.parameters['v'], self.parameters['deltav'], self.parameters['beta'], self.parameters['a'],pp=pre))
            
            if (error):            
                self.logger.info("{0:{pp}}  {1:{pp}}  {2:{pp}}  {3:{pp}}  {4:{pp}}  {5:{pp}}  {6:{pp}}  {7:{pp}}".format(self.error['B'], self.error['thB'], self.error['phiB'], self.error['tau'], 
                self.error['v'], self.error['deltav'], self.error['beta'], self.error['a'],pp=pre))


    def synthazel(self,method,stokes=None, returnRF=False, nlte=None,epsout=None, etaout=None, stimout=None):
        """
        Carry out the synthesis and returns the Stokes parameters directly from python user main program.
        ----------Parameters----------
        method = synthesis method for Hazel
        stokes : float
        An array of size [4 x nLambda] with the input Stokes parameter.                
        -------Returns-------
        stokes : float
        Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (4,nLambda)        
        """

        if (self.working_mode == 'inversion'):
            self.nodes_to_model()            
            self.to_physical()
        
        self.spectrum.synmethods.append(method) #here method is already a number

        B1In = np.asarray([self.parameters['B'],self.parameters['thB'], self.parameters['phiB']])
        hIn = self.height
        tau1In = self.parameters['tau']
        anglesIn = self.spectrum.los
        transIn = self.line_to_index[self.active_line]  #This was defined in en add_spectral/um
        mltp=self.spectrum.multiplets #only for making code shorter below
        lambdaAxisIn = self.wvl_axis - mltp[self.active_line]        
        nLambdaIn = len(lambdaAxisIn)
                        
        '''
        # Renormalize nbar so that its CLV is the same as that of Allen, but with a decreased I0
        # If I don't do that, fitting profiles in the umbra is not possible. The lines become in
        # emission because the value of the source function, a consequence of the pumping radiation,
        # is too large. In this case, one needs to use beta to reduce the value of the source function.
        ratio = boundaryIn[0,0] / i0_allen(mltp[self.active_line], self.spectrum.mu)

        '''
        '''-----------EDGAR: THOUGHTS CONCERNING BOUNDARY CONDITIONS------------------------------
        
        Concerning the above comment and the following code...
        1) In the first layer of the transfer (where boundary cond. are applied), the value of ratio 
        ends up being the value of I set in boundary=[1,0,0,0], because one first multiplies the intensity
        boundary component by i0Allen, and later divides it again in ratio, hence it has no effect other
        than setting ratio to boundary[0]=1. This 1 is a normalization corresponding to the value 
        that the user wants to assume for the continuum intensity(e.g. 1 for Allen continuum at a given mu, 
        or smaller values for the lower continuum intensities such as that of an umbra). 
        Then, this value should be inferred or even inverted for fitting observations, 
        but in synthesis we provide it.

        2) The next block of code also sets the value of BoundaryIn for Hazel. The value of BoundaryIn 
        entering Hazel in the very first boundary layer was here the i0Allen with spectral dependence 
        because boundary=[1,0,0,0] is multiplied by i0Allen and broadcasted to wavelength in the spectrum
        methods subroutines. In upper layers the stokes in boundaryIn is the result of previous transfer 
        but the intensity value at the wings (boundaryIn[0,0,]) used to define 'ratio' is still the continuum intensity
        set to 1 by the user in the boundary keyword because transfer does not affect the far wings in absence
        of continuum opacity.

        3) Then, one could think of avoiding to multiply and divide by i0Allen by just setting ratio 
        directly to the spectrum.boundary for intensity set by the user(which is the fraction, normally 1, 
        of the i0Allen at mu; in a dark background, of course this number should be lower, but this is set 
        ad hoc by the user or ideally by the inversion code.)
        However, we have to multiply and divide by I0Allen in every step in order to use the updated boundaryIn of the 
        layer. Improving this is not important given the approximations associated to the anisotropy (see below).

        4)We want the possibility of defining the boundaryIn for Hazel as the spectrum of i0Allen 
        or, more generally, as spectral profiles introduced by the user for every Stokes parameter 
        (for instance coming from observations or from ad hoc profiles for experiments). 

        So we reuse the keyword boundary to allow the possibility 
        of containing a provided spectral dependence or just 4 numbers (that then would be broadcasted
        in wavelength as constant profiles as before).
        Complementing the boundary keyword, we define the new keyword i0fraction,
        which is just one number representing the fraction of i0Allen continuum for intensity 
        and for normalization of Stokes vector to Icontinuum, but this keyword .

        It seems correct to modify nbar as of the number of photons of continuum relative to I0Allen,
        but this is just an approximation limited by how Hazel includes the anisotropy. 
        The nbar does not belong to a layer, but to the rest of the 
        atmopshere that illuminates the layer because nbar is the number of photons used to estimate radiation field tensors in
        the SEE. One then can argue that nbar for only one layer is related to the background continuum
        but, for multiple layers, the lower layers can modify the nbar of subsequent upper layers, which 
        is inconsistent with using the same Allen anisotropy for all layers. 
        Hence we would like to estimate this variation, realizing also that the nbar is not the same 
        for line center than for the wings. As this theory is limited to zero order we just 
        use the continuum because the spectral variation is not considered (flat spectrum aproximation).
        The above inconsistency is an intrinsic limitation of Hazel that is directly associated to the way 
        we introduce the radiation field tensor (the anisotropy) in the calculations. In more general cases 
        we calculate the anisotropy layer by layer using the intensity(and polarization) that arrives to the layer after performing the 
        radiative transfer along surrounding shells, and thus the nbar comes implicitly set by the sourroundings 
        and modulated in the trasnfer.  
           
        Here, ratio only affects nbar (hence anisotropies), but not the transfer of stokes. 
        For the first layer, ratio=boundaryIn[0,0]/i0_allen gives a fraction without units as required for 
        ratio and nbar.boundaryIn[0,0] is later updated from previous output and the code does that fresh 
        division again for every subsequent layer. As boundaryIn[0,0] is the continuum intensity, it does not change 
        appreciably with the transfer in absence of continuum opacity (Hazel still does not have it anyways).
        But if the continuum opacity is introduced or the opacity of the previous layer
        was large in the wings, boundaryIn[0,0] decreases along the LOS, which decreases ratio too as the transfer advances.
        The only problem with this is that the reduction of ratio (hence of nbar) is the same for all 
        rays in the radiation field sphere, so anisotropic transfer is not considered,as explained above.
        In any case, as boundaryIn always has physical units in every step of the transfer, it is reasonable 
        to divide by I0Allen to get the reamaining number of photons per mode (the fraction nbar) at every layer. 
        
        '''
        #-------------------------------------------------
        #we multiply by i0Allen to get units right. When introducing ad hoc the boundary 
        #with spectral dependence, we shall do it normalized to I0Allen(instead of with physical units), 
        #so still multiplication by I0Allen is necessary here.
        #self.spectrum.boundary arrives already multiplied by i0fraction if necessary. 
        if (stokes is None):
            boundaryIn  = np.asfortranarray(np.zeros((4,nLambdaIn)))
            boundaryIn[0,:] = i0_allen(mltp[self.active_line], self.spectrum.mu) #hsra_continuum(mltp[self.active_line]) 
            boundaryIn *= self.spectrum.boundary[:,self.wvl_range[0]:self.wvl_range[1]]
        else:            
            boundaryIn = np.asfortranarray(stokes)

        '''
        EDGAR: the value of boundary that enters here 
        must be Ibackground(physical units)/I0Allen (all input and output Stokes shall always be
        relative to the Allen I0 continuum). If the first value of this quantity is 1.0 then we have 
        an Allen background. Otherwise, that first value represents the true background 
        continuum realtive to I0Allen which can be be a fraction  of I0Allen (i.e., is the i0fraction
        introduced in model.py as tentative keyword). 

        For spectrum.boundary=1, the boundary intensity given to Hazel synthesize routine should be
        Iboundary = 1 *I0Allen(lambda), such that it has physical units (is not relative).

        spectrum.boundary is I0(physical)/I0Allen
        boundaryIn is then spectrum.boundary*I0Allen = I0(physical)
        then ratio=boundaryIn/I0Allen=I0(physical)/I0Allen , 
        which is a fraction of 1.0 (relative to the I0llen), as desired for nbarIn.
        '''
        ratio = boundaryIn[0,0]/ i0_allen(mltp[self.active_line], self.spectrum.mu)

        #nbarIn are reduction factors of nbar Allen for every transition! This means that it has 4 elements
        #for Helium and 2 for Sodium for instance, but this number was hardcoded to 4.
        #In addition omegaIn was wrong because it was put to zero, meaning no anisotropy,
        #while ones would mean that we use Allen for these pars.
        #when different than 0.0 and 1.0 they are used as modulatory factors in  hazel
        nbarIn = self.nbar * ratio #np.ones(self.ntr) * ratio
        omegaIn = self.j20f #np.ones(self.ntr)    #Not anymore np.zeros(4) 
        j10In = self.j10   #remind j10 and j20f are vectors (one val per transition).

        betaIn = self.parameters['beta']      

        #-------------------------------------------------

        dopplerWidthIn = self.parameters['deltav']
        dampingIn = self.parameters['a']
        dopplerVelocityIn = self.parameters['v']

        #Check where self.index is updated. It is index of current chromosphere,from 1 to n_chromospheres. 
        args = (self.index, method, B1In, hIn, tau1In, boundaryIn, transIn, anglesIn, nLambdaIn,
            lambdaAxisIn, dopplerWidthIn, dampingIn, j10In, dopplerVelocityIn,
            betaIn, nbarIn, omegaIn, self.atompol,self.magopt,self.stimem,self.nocoh,np.asarray(self.dcol) )
        
        #2D opt coeffs yet (not height dependent), for current slab self.index
        l,stokes,epsout,etaout,stimout,error = hazel_code._synth(*args)

        if (error == 1):raise NumericalErrorHazel()

        ff = self.parameters['ff'] #include it in the return below
        
        return ff * stokes, epsout,etaout,stimout,error #/ hsra_continuum(mltp[self.active_line])