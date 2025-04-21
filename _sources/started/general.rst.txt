.. include:: ../hazel_name

General description
===================

The main working horse in |hazel2| is the ``Model`` class. A ``Model`` defines all
the necessary information used by |hazel2| and it is the only class that needs to be
used for any synthesis or inversion. Although it gives an enormous flexibility for fine-tuning
it makes it necessary to code a few lines if you want to do a synthesis or inversion. However,
do not panic, we'll explain all the steps in this documentation. The first step is to import
the package and create the class:

::

    import hazel
    mod = hazel.Model(config=None, working_mode='synthesis', verbose=0, randomization=None, root='')

The class ``Model`` admits the following keywords:

* ``config`` (optional): defines the configuration file to be used. The structure of this file is described in :ref:`configuration` and it is the standard way of usage, specially when doing inversions. If ``None`` or absent, no configuration file will be used.
* ``working_model`` : ``synthesis`` for synthesis mode and ``inversion`` for inversion mode.
* ``verbose`` (optional, default is 0) : verbosity level. No verbosity is 0. Increasing levels of verbosity are 1-4, each one giving more information, useful for debugging. Use 0 when doing large scale inversions.
* ``randomization`` (optional, default is None) : when carrying out inversions, it is useful to carry out several optimizations starting from different regions of the space of parameters and check for convergence. This defines the number of randomizations used.
* ``root`` (optional, default is '') : set the root directory where all models used by the configuration file will be placed.


Once the ``Model`` class has been intantiated, one can interact with it.