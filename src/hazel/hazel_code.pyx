# cython: language_level=3
from numpy cimport ndarray as ar
from numpy import empty, linspace, zeros, array

cdef extern:
	void c_hazel(int* index, int* synMethInput, double* B1Input, double* hInput, double* tau1Input, 
		double* boundaryInput, int* transInput, double* anglesInput, int* nLambdaInput, double* lambdaAxisInput,
		double* dopplerWidthInput, double* dampingInput, double* j10Input, double* dopplerVelocityInput, 
		double* betaInput, double* nbarInput, double* omegaInput, 
		int* atompolInput,int* magoptInput,int* stimemInput,int* nocohInput, double* dcolInput,
		double* wavelengthOut, double* stokesOut,double* epsOut,double* etaOut,double* stimOut, int* error)
		
	void c_init(int* nchar,char* atomfile, int* verbose,int* ntransOutput) #EDGAR: added nchar, atomfile, and output par ntrans
	void c_exit(int* index)

#routine called by python synthazel
def _synth(int index=1, int synMethInput=5, ar[double,ndim=1] B1Input=zeros(3), double hInput=3.0, 
	double tau1Input=1.0, 
	ar[double,ndim=2,mode='fortran'] boundaryInput=zeros((4,128)), int transInput=1, ar[double,ndim=1] anglesInput=zeros(3), 
	int nLambdaInput=128, ar[double,ndim=1] lambdaAxisInput=linspace(-1.5,2.5,128),  
	double dopplerWidthInput=5.0, double dampingInput=0.0, ar[double,ndim=1] j10Input=zeros(4),double dopplerVelocityInput=0.0, 
	double betaInput=1.0, ar[double,ndim=1] nbarInput=zeros(4), 
	ar[double,ndim=1] omegaInput=zeros(4),
	int atompolInput=1,int magoptInput=1,int stimemInput=1,int nocohInput=0, 
	ar[double,ndim=1] dcolInput=zeros(3)):
	

	#EDGAR: 
	#nLambdaInput, lambdaAxisInput,nLambdaInput are hardcoded to 128 points or is just an initialization?
	#same for nbarInput and omegaInput which seems initizialized to 4 transitions.
	"""
	Carry out a synthesis with Hazel
	
	Args: (see the manual for the meaning of all of them)
		index: (int) index of atmosphere
		B1Input: (float) vector of size 3 with the magnetic field vector in spherical coordinates for the first component
		hInput: (float) height
		tau1Input: (float) optical depth of the first component
		boundaryInput: (float) vector of size 4xnLambda with the boundary condition for (I,Q,U,V)
		transInput: (int) transition to compute from the model atom
		anglesInput: (float) vector of size 3 describing the LOS
		lambdaAxisInput: (float) vector of size 2 defining the left and right limits of the wavelength axis
		nLambdaInput: (int) number of wavelength points
		dopplerWidth1Input: (float) Doppler width of the first component
		dampingInput: (float) damping
		dopplerVelocityInput: (float) bulk velocity affecting the first component
		betaInput: (float) enhancement factor for the source function of component 1 to allow for emission lines in the disk
		nbarInput: (float) vector of size 4 to define nbar for every transition of the model atom (set them to zero to use Allen's)
		omegaInput: (float) vector of size 4 to define omega for every transition of the model atom (set them to zero to use Allen's)
		
    Returns:
        wavelengthOutput: (float) vector of size nLambdaInput with the wavelength axis
        stokesOutput: (float) array of size (4,nLambdaInput) with the emergent Stokes profiles
        epsOutput: (float) array of size (4,nLambdaInput) with the emissivity vector at each wavelength
        etaOutput: (float) array of size (7,nLambdaInput) with the independent elements of K matrix at each wavelength
		error: (int) zero if everything went OK
	"""
	
	cdef:
		ar[double,ndim=1] wavelengthOut = empty(nLambdaInput, order='F')
		ar[double,ndim=2] stokesOut = empty((4,nLambdaInput), order='F')
		ar[double,ndim=2] epsOut = empty((4,nLambdaInput), order='F')
		ar[double,ndim=2] etaOut = empty((7,nLambdaInput), order='F')
		ar[double,ndim=2] stimOut = empty((7,nLambdaInput), order='F')
		int error

	#calls fortran routine c_hazel in hazel_py.f90
	c_hazel(&index, &synMethInput, &B1Input[0], &hInput, &tau1Input,  
		&boundaryInput[0,0], &transInput, &anglesInput[0], &nLambdaInput, &lambdaAxisInput[0],  
		&dopplerWidthInput, &dampingInput, &j10Input[0], &dopplerVelocityInput, 
		&betaInput, &nbarInput[0], &omegaInput[0], 
		&atompolInput,&magoptInput,&stimemInput,&nocohInput,&dcolInput[0],
		<double*> wavelengthOut.data, <double*> stokesOut.data,<double*> epsOut.data, <double*> etaOut.data,
		<double*> stimOut.data, &error)
    
	return wavelengthOut, stokesOut, epsOut, etaOut, stimOut, error
	
def _init(str atomfile, int verbose=0):
	"""
	Initialize and do some precomputations that can be avoided in the subsequent calls to the synthesis
	Args:
        atomfile: (str) name of the input atom file .mod to be read. This is a C string
    Returns:
        None

	EDGAR:The line atomfile.encode() converts the input atomfile string to utf8

	"""
	ftmp = atomfile.encode()
	cdef:
		int nchar = len(atomfile)
		char* atomfileInput = ftmp
		#int ntransOutput = 0
		int ntransOutput = 0

	
	#c_init(&nchar, &atomfileInput[0], &verbose) #&ntransOutput
	#this calls init routine in hazel_py.f90 and return ntrans 
	c_init(&nchar, &atomfileInput[0], &verbose, &ntransOutput)

	return ntransOutput #EDGAR:now it returns ntrans when called from model.py

def _exit(int index):
		
	c_exit(&index)

