.. include:: ../hazel_name

Generators
==========

Although you can follow what you've learnt in the previous section to generate
the input files, |hazel2| can generate the input files for you in the appropriate
format.

Observations
------------

You can generate the appropriate files with the observations in 1D or 3D formats.
Once you have the data in `numpy` arrays, you can instantiate a `File_observation`
as follows for a 1D observation:

::

    tmp = hazel.tools.File_observation(mode='single')
    tmp.set_size(n_lambda=128, n_pixel=1)

or as follows for a 3D observation:

::
 
    tmp = hazel.tools.File_observation(mode='multi')
    tmp.set_size(n_lambda=128, n_pixel=10) 

In this case, `tmp` contains a dictionary `obs` with the following keys with the correct dimensions:

::

    In [1]: tmp.obs.keys()
    Out[1]: dict_keys(['stokes', 'sigma', 'los', 'boundary', 'wavelength', 'weights', 'mask'])

You can then fill up these arrays with the appropriate data and finally save the observations
with:

::

    tmp.save('test')

This will generate the files with the Stokes profiles, the wavelength and the wavelength weights, all
of them with the same name and different extensions. These files can then be used with |hazel2| for
carrying out inversions.

Photospheres
------------

Something similar can be achieved for photospheric models. First instantiate the following class:

::

    tmp = hazel.tools.File_photosphere(mode='single')

The object `tmp` contains now a dictionary with the following keys:

::

    In [1]: tmp.model.keys()
    Out[1]: dict_keys(['model', 'ff'])

You can then fill them appropriately with the depth stratification of the model or you can
use the `set_default` method of the object to fill it up with standard model atmospheres. For
the moment, you can just use the HSRA model, but more models will come in the future:

::

    tmp.set_default(n_pixel=1, default='hsra')
    tmp.save('photosphere')

To get a list of possible models, just type:

::

    tmp.list_models()


In the multipixel case, just use:

::

    tmp = hazel.tools.File_photosphere(mode='multi')
    tmp.set_default(n_pixel=10, default='hsra')
    tmp.save('photosphere')

Chromospheres
-------------

The same applies to chromospheres. An example follows for a multipixel atmosphere:

::

    tmp = hazel.tools.File_chromosphere(mode='multi')
    tmp.set_default(n_pixel=10, default='disk')
    tmp.save('chromosphere')

The options for the defaults are `disk` for on-disk observations and `offlimb` for off-limb observations.
