import numpy as np
from scipy.io import netcdf

def genNetCDF(wavelength, stI, stQ, stU, stV, sigmaI, sigmaQ, sigmaU, sigmaV, boundary, height, obsTheta, obsGamma, mask, pars, normalization, outputFile):
    """
    This routine generates a NetCDF file with the observations ready for Hazel-MPI
    
    Args:
        wavelength (float): array of size [nlambda]
        stI (float): array of size [npixel, nlambda] with Stokes I
        stQ (float): array of size [npixel, nlambda] with Stokes Q
        stU (float): array of size [npixel, nlambda] with Stokes U
        stV (float): array of size [npixel, nlambda] with Stokes V
        sigmaI (float): array of size [npixel, nlambda] with the noise in Stokes I
        sigmaQ (float): array of size [npixel, nlambda] with the noise in Stokes Q
        sigmaU (float): array of size [npixel, nlambda] with the noise in Stokes U
        sigmaV (float): array of size [npixel, nlambda] with the noise in Stokes V
        boundary (float): array of size [npixel, 4] with the boundary conditions [I0,Q0,U0,V0] for every pixel
        height (float): array of size [npixel] indicating the height of the pixel over the surface in arcsec
        obsTheta (float): array of size [npixel] indicating the angle of the observation in degrees 
        obsGamma (float): array of size [npixel] the angle of the reference for Stokes Q
        mask (float): array of the original dimensions of the observations that is used later to reconstruct the inverted maps [nx,ny]
        pars (float): array of size [nparameters,npixel] that gives the initial value of the parameters
                            The size depends on the radiative transfer option:
                                * 1-component (vector of size 8): B, thetaB, chiB, tau, vdop, a, vmac, beta
                                * 2-component 1+1 with same field (vector of size 11): B, thetaB, chiB, tau1, tau2, vdop, a, vmac1, vmac2, beta, beta2
                                * 2-component 1+1 with different field (vector of size 15): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, beta, beta2
                                * 2-component 2 with different field with ff (vector of size 16): B1, thetaB1, chiB1, B2, thetaB2, chiB2, tau1, tau2, vdop1, vdop2, a, vmac1, vmac2, ff, beta, beta2
        normalization (string): 'continuum' or 'peak', indicating how the profiles are normalized
        outputFile (float): output file
    """
    
    nPixel, nLambda = stI.shape
    nCols = 8
    
    obsMap = np.zeros((8,nLambda,nPixel))
    
    obsMap[0,:,:] = stI.T
    obsMap[1,:,:] = stQ.T
    obsMap[2,:,:] = stU.T
    obsMap[3,:,:] = stV.T
    obsMap[4,:,:] = sigmaI.T
    obsMap[5,:,:] = sigmaQ.T
    obsMap[6,:,:] = sigmaU.T
    obsMap[7,:,:] = sigmaV.T
    
    obsMap = np.transpose(obsMap,axes=(2,1,0))
    
    dimMap = mask.shape
                
# Variable dimensions
    fileID = netcdf.netcdf_file(outputFile, 'w')
    nPixDim = fileID.createDimension('npixel', nPixel)
    nColDim = fileID.createDimension('ncolumns', nCols)
    nStokesParDim = fileID.createDimension('nstokes_par', 4)
    nParsDim = fileID.createDimension('nparameters', pars.shape[0])
    nLambdaDim = fileID.createDimension('nlambda', nLambda)
    nXDim = fileID.createDimension('nx', dimMap[0])
    nYDim = fileID.createDimension('ny', dimMap[1])

# Variable definition
# Remember that variables are written in C format, so that they are reversed with respect to Fortran
    lambdaID = fileID.createVariable('lambda','f8',('nlambda',))
    stokesID = fileID.createVariable('map','f8',('npixel','nlambda','ncolumns',))
    boundaryID = fileID.createVariable('boundary','f8',('npixel','nstokes_par',))
    heightID = fileID.createVariable('height','f8',('npixel',))
    obsThetaID = fileID.createVariable('obs_theta','f8',('npixel',))
    obsGammaID = fileID.createVariable('obs_gamma','f8',('npixel',))
    maskID = fileID.createVariable('mask','i2',('nx','ny',))
    parsInitID = fileID.createVariable('pars_initial','f8',('npixel','nparameters',))
    normalizationID = fileID.createVariable('normalization','f8',())


    lambdaID[:] = wavelength
    stokesID[:] = obsMap
    boundaryID[:] = boundary
    heightID[:] = height
    obsThetaID[:] = obsTheta
    obsGammaID[:] = obsGamma
    maskID[:] = mask
    parsInitID[:] = pars.T
    if (normalization == 'peak'):
        normalizationID.assignValue(1)
    if (normalization == 'continuum'):
        normalizationID.assignValue(0)
    
    fileID.close()