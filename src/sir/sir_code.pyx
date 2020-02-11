# cython: language_level=3
from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray
import numpy as np
cimport numpy as np

cdef extern:
	void c_init_externalfile(int *index, int *nchar, char *file_lines, int *nLambda)
	void c_init(int *index, int *nblend, int *lines, int *atom, int *istage, double *wlength, double *zeff, double *energy, double *loggf,
		int *mult1, int *mult2, int *design1, int *design2, double *tam1, double *tam2, double *alfa, double *sigma,
		double *lambda0, double *lambda1, int *nsteps)
	void c_setpsf(int *nPSF, float *xPSF, float *yPSF)
	void c_synthrf(int *index, int *nDepth, int *nLambda, double *macroturbulence, double *model, double *stokes, double *rt, double *rp, double *rh,
		double *rv, double *rf, double *rg, double *rm, double *rmac, int *error)
	void c_synth(int *index, int *nDepth, int *nLambda, double *macroturbulence, double *model, double *stokes, int *error)

def init_externalfile(int index, str file):
	cdef:
		int nLambda
		ftmp = file.encode('UTF-8')
		cdef char* file_lines = ftmp
		cdef int   nchar      = len(file)
		
	c_init_externalfile(&index, &nchar, file_lines, &nLambda)

	return nLambda

def init(int index, int nblend, ar[int, ndim=1] lines, ar[int, ndim=1] atom, ar[int, ndim=1] istage, ar[double, ndim=1] wlength, ar[double, ndim=1] zeff, 
	ar[double, ndim=1] energy, ar[double, ndim=1] loggf, ar[int, ndim=1] mult1, ar[int, ndim=1] mult2, 
	ar[int, ndim=1] design1, ar[int, ndim=1] design2, ar[double, ndim=1] tam1, ar[double, ndim=1] tam2, 
	ar[double, ndim=1] alfa, ar[double, ndim=1] sigma, double lambda0_in, double lambda1_in, int n_steps_in):
	
	
	c_init(&index, &nblend, &lines[0], &atom[0], &istage[0], &wlength[0], &zeff[0], &energy[0], &loggf[0], &mult1[0], &mult2[0],
		&design1[0], &design2[0], &tam1[0], &tam2[0], &alfa[0], &sigma[0], &lambda0_in, &lambda1_in, &n_steps_in)

	return 

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
		int error
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

	c_synth(&index, &nDepth, &nLambda, &macroturbulence, &model[0,0], <double*> stokes.data, &error)
	
	return stokes, error

def synthRF(int index, int nLambda, ar[double, ndim=1] log_tau, ar[double, ndim=1] T, ar[double, ndim=1] Pe, ar[double, ndim=1] vmic, 
	ar[double, ndim=1] vlos, ar[double, ndim=1] Bx, ar[double, ndim=1] By, ar[double, ndim=1] Bz, 
	double macroturbulence):

	cdef:
		int nDepth = len(log_tau)
		int error
		ar[double, ndim=2] model = empty((8,nDepth), order='F')
		ar[double, ndim=2] stokes = empty((5,nLambda), order='F')
		ar[double, ndim=3] rt = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rp = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rh = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rv = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rf = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rg = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=3] rm = empty((4,nLambda,nDepth), order='F')
		ar[double, ndim=2] rmac = empty((4,nLambda), order='F')
		
	model[0,:] = log_tau
	model[1,:] = T
	model[2,:] = Pe
	model[3,:] = vmic
	model[4,:] = np.sqrt(Bx**2 + By**2 + Bz**2)
	model[5,:] = vlos
	model[6,:] = 180.0 / np.pi * np.arccos(Bz / (model[4,:] + 1e-8))
	model[7,:] = 180.0 / np.pi * np.arctan2(By, Bx)

	c_synthrf(&index, &nDepth, &nLambda, &macroturbulence, &model[0,0], <double*> stokes.data, <double*> rt.data, <double*> rp.data, 
		<double*> rh.data, <double*> rv.data, <double*> rf.data, <double*> rg.data, <double*> rm.data, 
		<double*> rmac.data, &error)
	
	return stokes, [rt, rp, rh, rv, rf, rg, rm, rmac], error