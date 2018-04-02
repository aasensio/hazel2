from numpy cimport ndarray as ar
from numpy import empty, linspace, zeros, array

cdef extern:
	void c_hazel(int* index, double* B1Input, double* hInput, double* tau1Input, 
		double* boundaryInput, int* transInput, double* anglesInput, int* nLambdaInput, double* lambdaAxisInput,
		double* dopplerWidthInput, double* dampingInput, double* dopplerVelocityInput, 
		double* betaInput, double* nbarInput, double* omegaInput, 
		double* wavelengthOutput, double* stokesOutput)
		
	void c_init()

def _synth(int index=1, ar[double,ndim=1] B1Input=zeros(3), double hInput=3.0, 
	double tau1Input=1.0, 
	ar[double,ndim=2,mode='fortran'] boundaryInput=zeros((4,128)), int transInput=1, ar[double,ndim=1] anglesInput=zeros(3), 
	int nLambdaInput=128, ar[double,ndim=1] lambdaAxisInput=linspace(-1.5,2.5,128),  
	double dopplerWidthInput=5.0, double dampingInput=0.0, double dopplerVelocityInput=0.0, 
	double betaInput=1.0, ar[double,ndim=1] nbarInput=zeros(4), 
	ar[double,ndim=1] omegaInput=zeros(4)):
	
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
	"""
	
	cdef:
		ar[double,ndim=1] wavelengthOutput = empty(nLambdaInput, order='F')
		ar[double,ndim=2] stokesOutput = empty((4,nLambdaInput), order='F')
		
	c_hazel(&index, &B1Input[0], &hInput, &tau1Input,  
		&boundaryInput[0,0], &transInput, &anglesInput[0], &nLambdaInput, &lambdaAxisInput[0],  
		&dopplerWidthInput, &dampingInput, &dopplerVelocityInput, 
		&betaInput, &nbarInput[0], &omegaInput[0], <double*> wavelengthOutput.data, 
		<double*> stokesOutput.data)
    
	return wavelengthOutput, stokesOutput
	
def _init():
	"""
	Initialize and do some precomputations that can be avoided in the subsequent calls to the synthesis
	
	Args:
        None
    Returns:
        None
	"""
	c_init()
