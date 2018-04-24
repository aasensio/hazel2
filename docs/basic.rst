.. include:: hazel_name
Basic usage
===========

Programmatically
----------------

For simple calculations, like synthesizing spectral lines in simple models,
|hazel2| can be used in programmatic mode. For instance, let us generate a spectral
window in the near-infrared and synthesize the He I 10830 A line with some
parameters.

::
    
    import hazel
    import matplotlib.pyplot as plt
    mod = hazel.Model(working_mode='synthesis')
    mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, 150], 'topology': 'ch1',
	'LOS': [0.0,0.0,90.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})
    mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 3.0, 'Line': '10830', 'Wavelength': [10826, 10833]})
    mod.setup()

    f, ax = plt.subplots(nrows=2, ncols=2)
    ax = ax.flatten()

    for j in range(5):
        mod.atmospheres['ch1'].set_parameters([0.0,0.0,100.0*j,1.0,0.0,8.0,1.0,0.0],1.0)
        mod.synthesize()

        for i in range(4):
            ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
    plt.show()
    

As you see, we first generate the Hazel model for synthesis. We then create a spectral window with parameters
passed as a dictionary. Then we add a chromosphere and finally call `setup` to finalize the model setup.
We then change the model parameters of the chromosphere and synthesize the spectrum.
Access to the synthesized spectra is done on the ``stokes`` key of the ``mod.spectrum`` dictionary.
All the details of the dictionaries to pass and how to generate more complicated
atmospheres can be found in :ref:`programmatically`.

With configuration file
-----------------------

Perhaps the easiest way of running |hazel2| is through the human-friendly configuration
files described in :ref:`configuration`. The configuration files used in this section
are present in `<https://github.com/aasensio/hazel2/examples>`_

Single pixel mode
^^^^^^^^^^^^^^^^^

Calculations in single-pixel mode (only one pixel synthesis/inversion) are very easy
to do. The following code uses a configuration file with `1d` inputs files, 
whose format is described in :ref:`input`. The following one
carries out synthesis:

::

    mod = hazel.Model('conf_single_syn.ini')
    mod.open_output()
    mod.synthesize()
    mod.write_output()
    mod.close_output()

and this one carries out the inversion:

::

    mod = hazel.Model('conf_single_inv.ini')
    mod.read_observation()
    mod.open_output()
    mod.invert()
    mod.write_output()
    mod.close_output()

Summarizing, just create the model using the configuration file and then 
call synthesize/invert. Note that we open and close the output because
we want to generate a file with the synthetic profiles or the model
parameters of the inversion.

Serial mode
^^^^^^^^^^^

When many pixels need to be synthesized/inverted, one can pass appropriate input
multidimensional files that can deal with many pixels. They can be synthesized
in serial mode (using only one CPU) with the following code. It makes use of
an `iterator`, an object that  carries out the work for all the pixels in the
input files.

::

    iterator = hazel.iterator(use_mpi=False)    
    mod = hazel.Model('conf_nonmpi_syn1d.ini')
    iterator.use_model(model=mod)
    iterator.run_all_pixels()

This case is slightly more complicated. We need to instantiate an iterator, telling it
not to use MPI. The, we instantiate the model and pass it to the iterator, which is
then used to run through all pixels.

MPI mode
^^^^^^^^

When working with many pixels, many-core computers can be used. This will use a
master-slave approach, in which the master sends pixels to each one of the available
slaves to do the work. In this case, we tell the iterator to use MPI and act differently 
for the master (rank=0) or the salves (rank>0). Note that only the master reads the
configuration file, and it will be broadcasted to all slaves internally.

::

    iterator = hazel.iterator(use_mpi=True)
    rank = iterator.get_rank()

    if (rank == 0):    
        mod = hazel.Model('conf_mpi_invh5.ini')
        iterator.use_model(model=mod)
    else:
        iterator.use_model()

    iterator.run_all_pixels()
