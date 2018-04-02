Installation
============

Hazel v2.0 is a Python module with some wrapped compiled Fortran routines. It should
be pretty easy to install once you have a working compiler in your system.

There are different ways to install Hazel v2.0, but the best is to install
it into a virtual environment either with `pip <http://www.virtualenv.org/>`_ or `conda <https://conda.io/docs/user-guide/tasks/manage-environments.html/>`_.

From source
-----------
It can also be installed from sources by cloning this repository and installing it:

::
    
    git clone https://github.com/aasensio/hazel
    python setup.py install

Remember that if you want to develop the code, it is useful to use instead:

::
    
    python setup.py develop

Using pip
---------

Hazel v2.0 can be easily installed from PyPI by invoking:

::

    pip install hazel

Using conda
-----------

It can also be installed from the Anaconda Cloud by invoking:

::

    conda install hazel

Requirements
------------
Hazel v2.0 depends on the following external packages, that should be
pretty straightforward to install:

* numpy
* h5py
* scipy
* astropy
* mpi4py
* configobh

All of them can be installed in Anaconda with:

::

    conda install numpy h5py scipy astropy mpi4py configobj