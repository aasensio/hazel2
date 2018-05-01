.. include:: hazel_name
Installation
============

|hazel2| is a Python 3 module with some wrapped compiled Fortran routines. It should
be pretty easy to install once you have a working compiler in your system.

There are different ways to install |hazel2|, but the best is to install
it into a virtual environment either with `pip <http://www.virtualenv.org/>`_ or `conda <https://conda.io/docs/user-guide/tasks/manage-environments.html/>`_,
which makes everything much more safer, plus making sure that all packages are installed for the code.
For example, once you have installed `Miniconda <https://conda.io/miniconda.html/>`_, you can generate
a new environment and install the dependencies (you can install whatever version of Python 3 you desire):

::

    conda create -n hazel_env python=3.6 numpy h5py scipy astropy mpi4py configobj tqdm
    source activate hazel_env

If you also want to use the GUI, you need to add two new dependencies:

::

    conda create -n hazel_env python=3.6 
    source activate hazel_env
    conda install -c conda-forge cython numpy scipy h5py sphinx sphinx_rtd_theme tqdm astropy configobj matplotlib numpydoc mpich mpi4py    



From source
-----------
It can also be installed from sources by cloning this repository and installing it:

::
    
    git clone https://github.com/aasensio/hazel2
    python setup.py install

Remember that if you want to develop the code, it is useful to use instead:

::
    
    python setup.py develop

and when you want to uninstall the development version, type:

::
    
    python setup.py develop --uninstall

..
    Using pip
    ---------

    |hazel2| can be easily installed from PyPI by invoking:

    ::

        pip install hazel

    Using conda
    -----------

    It can also be installed from the Anaconda Cloud by invoking:

    ::

        conda install hazel

Requirements
------------
|hazel2| depends on the following external packages, that should be
pretty straightforward to install:

* ``numpy``
* ``h5py``
* ``scipy``
* ``astropy``
* ``mpi4py``
* ``configobj``
* ``tqdm``

If you want to use the GUI, you also need to install the following dependencies:

* matplotlib
* pyqt5

All of them can be installed in Anaconda with:

::

    conda install numpy h5py scipy astropy mpi4py configobj
