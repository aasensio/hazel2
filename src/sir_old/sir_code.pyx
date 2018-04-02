from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray
import numpy as np

cdef extern:
	void c_init(int *index, int *nLambda)
	void c_setpsf(int *nPSF, float *xPSF, float *yPSF)
	void c_synthrf(int *index, int *nDepth, int *nLambda, float *macroturbulence, float *filling, float *stray, float *model, float *stokes, float *rt, float *rp, float *rh,
		float *rv, float *rf, float *rg, float *rm, float *rmac)
	void c_synth(int *index, int *nDepth, int *nLambda, float *macroturbulence, float *filling, float *stray, float *model, float *stokes)

def init(int index):
	cdef:
		int nLambda
	
	c_init(&index, &nLambda)

	return nLambda

def setPSF(ar[float, ndim=1] xPSF, ar[float, ndim=1] yPSF):
	cdef:
		int nLambdaPSF = xPSF.shape[0]

	c_setpsf(&nLambdaPSF, &xPSF[0], &yPSF[0])

	return

def synth(int index, int nLambda, ar[float, ndim=2] modelIn, float macroturbulence, float filling, float stray):

	cdef:
		int nDepth = modelIn.shape[0]
		ar[float, ndim=2, mode="c"] model		
		ar[float, ndim=2] stokes = empty((5,nLambda), order='F', dtype=np.float32)		
		
	# Make sure that the 2D array is C_CONTIGUOUS
	model = ascontiguousarray(modelIn)

	c_synth(&index, &nDepth, &nLambda, &macroturbulence, &filling, &stray, &model[0,0], <float*> stokes.data)
	
	return stokes

def synthRF(int index, int nLambda, ar[float, ndim=2] modelIn, float macroturbulence, float filling, float stray):

	cdef:
		int nDepth = modelIn.shape[0]
		ar[float, ndim=2, mode="c"] model		
		ar[float, ndim=2] stokes = empty((5,nLambda), order='F', dtype=np.float32)
		ar[float, ndim=3] rt = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rp = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rh = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rv = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rf = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rg = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rm = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=2] rmac = empty((4,nLambda), order='F', dtype=np.float32)
		
	# Make sure that the 2D array is C_CONTIGUOUS
	model = ascontiguousarray(modelIn)

	c_synthrf(&index, &nDepth, &nLambda, &macroturbulence, &filling, &stray, &model[0,0], <float*> stokes.data, <float*> rt.data, <float*> rp.data, 
		<float*> rh.data, <float*> rv.data, <float*> rf.data, <float*> rg.data, <float*> rm.data, <float*> rmac.data)
	
	return stokes, [rt, rp, rh, rv, rf, rg, rm, rmac]