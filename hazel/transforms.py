import numpy as np

__all__ = ['physical_to_transformed', 'transformed_to_physical']

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
    return np.log(x / (1.0 - x))

def inv_logit(self, x):
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
    return _logit( (x-lower) / (upper - lower))

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

def jacobianTransformedParameters(x):
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
    temp = self.inv_logit(x)
    return (upper - lower) * temp * (1.0 - temp)