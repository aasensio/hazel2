import numpy as np

__all__ = ['physical_to_transformed', 'transformed_to_physical']

_eps = 1e-8

def logit(x):
    """
    Logit function
    
    Parameters
    ----------
    x : float
        Any array

    Returns
    -------
    logit : float
        Logit transform of the input
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        out = np.log(x / (1.0 - x))
    return out

def inv_logit(x):
    """
    Inverse logit function
    
    Parameters
    ----------
    x : float
        Any array

    Returns
    -------
    inv_logit : float
        Inverse logit transform of the input
    """
    return 1.0 / (1.0 + np.exp(-x))

def physical_to_transformed(x, lower, upper):
    """
    Transform from physical parameters to unconstrained physical parameters
    
    Parameters
    ----------
    x : float
        Any array
    lower : float
        Lower limit of the parameter
    upper : float
        Upper limit of the parameter

    Returns
    -------
    out : float
        Transformed parameters
    """
    return logit( (x-lower) / (upper - lower))

def transformed_to_physical(x, lower, upper):
    """
    Transform from unconstrained physical parameters to physical parameters
    
    Parameters
    ----------
    x : float
        Any array
    lower : float
        Lower limit of the parameter
    upper : float
        Upper limit of the parameter

    Returns
    -------
    out : float
        Transformed parameters
    """
    return lower + (upper - lower) * inv_logit(x)

def jacobian_transformation(x, lower, upper):
    """
    Compute the Jacobian of the transformation from unconstrained parameters to physical parameters
    
    Parameters
    ----------
    x : float
        Any array
    lower : float
        Lower limit of the parameter
    upper : float
        Upper limit of the parameter

    Returns
    -------
    out : float
        Transformed parameters
    """
    temp = inv_logit(x)
    return (upper - lower) * temp * (1.0 - temp)

def jacobian_transformed_to_physical(x, lower, upper):
    """
    Compute the Jacobian of the transformation from unconstrained parameters to physical parameters
    
    Parameters
    ----------
    x : float
        Any array
    lower : float
        Lower limit of the parameter
    upper : float
        Upper limit of the parameter

    Returns
    -------
    out : float
        Transformed parameters
    """    
    return (upper - lower) / ((lower - x) * (x - upper))