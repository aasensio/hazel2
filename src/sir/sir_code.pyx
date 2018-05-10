from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray
import numpy as np

cdef extern:
	void c_init(int *index, int *nchar, char *file_lines, int *nLambda)
	void c_setpsf(int *nPSF, float *xPSF, float *yPSF)
	void c_synthrf(int *index, int *nDepth, int *nLambda, float *macroturbulence, float *model, float *stokes, float *rt, float *rp, float *rh,
		float *rv, float *rf, float *rg, float *rm, float *rmac)
	void c_synth(int *index, int *nDepth, int *nLambda, double *macroturbulence, double *model, double *stokes)

def init(int index, str file):
	cdef:
		int nLambda
		ftmp = file.encode('UTF-8')
		cdef char* file_lines = ftmp
		cdef int   nchar      = len(file)
		
	c_init(&index, &nchar, file_lines, &nLambda)

	return nLambda

def setPSF(ar[float, ndim=1] xPSF, ar[float, ndim=1] yPSF):
	cdef:
		int nLambdaPSF = xPSF.shape[0]

	c_setpsf(&nLambdaPSF, &xPSF[0], &yPSF[0])

	return

def synth(int index, int nLambda, ar[double, ndim=1] log_tau, ar[double, ndim=1] T, ar[double, ndim=1] Pe, ar[double, ndim=1] vmic, 
	ar[double, ndim=1] vlos, ar[double, ndim=1] Bx, ar[double, ndim=1] By, ar[double, ndim=1] Bz, 
	double macroturbulence):
	cdef:
		int nDepth = len(log_tau)
		ar[double, ndim=2] model = empty((8,nDepth), order='F')
		ar[double, ndim=2] stokes = empty((5,nLambda), order='F')
		
	model[0,:] = log_tau
	model[1,:] = T
	model[2,:] = Pe
	model[3,:] = vmic
	model[4,:] = np.sqrt(Bx**2 + By**2 + Bz**2)
	model[5,:] = vlos	
	model[6,:] = 180.0 / np.pi * np.arccos(Bz / (model[4,:] + 1e-8))      # Regularize in case B=0
	model[7,:] = 180.0 / np.pi * np.arctan2(By, Bx)

	c_synth(&index, &nDepth, &nLambda, &macroturbulence, &model[0,0], <double*> stokes.data)
	
	return stokes

def synthRF(int index, int nLambda, ar[float, ndim=1] log_tau, ar[float, ndim=1] T, ar[float, ndim=1] Pe, ar[float, ndim=1] vmic, 
	ar[float, ndim=1] vlos, ar[float, ndim=1] Bx, ar[float, ndim=1] By, ar[float, ndim=1] Bz, 
	float macroturbulence):

	cdef:
		int nDepth = len(log_tau)
		ar[float, ndim=2] model = empty((8,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=2] stokes = empty((5,nLambda), order='F', dtype=np.float32)
		ar[float, ndim=3] rt = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rp = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rh = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rv = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rf = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rg = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=3] rm = empty((4,nLambda,nDepth), order='F', dtype=np.float32)
		ar[float, ndim=2] rmac = empty((4,nLambda), order='F', dtype=np.float32)
		
	model[0,:] = log_tau
	model[1,:] = T
	model[2,:] = Pe
	model[3,:] = vmic
	model[4,:] = np.sqrt(Bx**2 + By**2 + Bz**2)
	model[5,:] = vlos
	model[6,:] = 180.0 / np.pi * np.arccos(Bz / model[4,:])
	model[7,:] = 180.0 / np.pi * np.arctan2(By, Bx)

	c_synthrf(&index, &nDepth, &nLambda, &macroturbulence, &model[0,0], <float*> stokes.data, <float*> rt.data, <float*> rp.data, 
		<float*> rh.data, <float*> rv.data, <float*> rf.data, <float*> rg.data, <float*> rm.data, <float*> rmac.data)
	
	return stokes, [rt, rp, rh, rv, rf, rg, rm, rmac]