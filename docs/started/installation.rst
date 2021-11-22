.. include:: ../hazel_name

Installation
============

Generalities
------------

|hazel2| is a Python 3 module with some wrapped compiled Fortran routines. It should
be pretty easy to install once you have a working compiler in your system.

There are different ways to install |hazel2|, but the best is to install
it into a virtual environment either with `pip <http://www.virtualenv.org>`_ or `conda <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_,
which makes everything much more safer, plus making sure that all packages are installed for the code.
For example, once you have installed `Miniconda <https://conda.io/miniconda.html>`_, you can generate
a new environment and install the dependencies (you can install whatever version of Python 3 you desire).
For consistency, it is better to use the gfortran compilers from Anaconda. In the
following we show how install all packages with Anaconda. Anyway, feel free to use
the gfortran compiler from your system but you might have some issues.

For a Linux OS, type:

::

    conda create -n hazel_env python=3.6
    conda activate hazel_env
    conda install -c conda-forge cython numpy h5py tqdm scipy astropy mpich mpi4py configobj gfortran_linux-64 gcc_linux-64 asciitree

For Mac OS, type:

::

    conda create -n hazel_env python=3.6
    conda activate hazel_env
    conda install -c conda-forge cython numpy h5py tqdm scipy astropy mpich mpi4py configobj gfortran_osx-64 gcc_osx-64 asciitree

The specific packages for the Fortran and C compilers might change depending on the version of Anaconda
that you have installed in your system. Their name might also be `gfortran` and `gcc`.
Feel free to add `matplotlib` for doing some plots. You will need it if you want to run some of the examples
shown in this documentation in your computer.

If you also want to use the GUI, you need to add two new dependencies:

::

    conda create -n hazel_env python=3.6 
    conda activate hazel_env
    conda install -c conda-forge cython numpy scipy h5py tqdm astropy configobj matplotlib mpich mpi4py pyqt asciitree

Remember to add `ipython` if you are using this shell to run |hazel2|. Otherwise, it will use
`ipython` from another environment and you can get confused.

If you want to use the NLTE neural option for synthesizing the Ca II 854.2 nm line, the packages ``pytorch``
and ``pytorch_geometric`` should be installed. Take a look at the documentation for `PyTorch <https://pytorch.org/>`_
and `PyTorch Geometric <https://pytorch-geometric.readthedocs.io/en/latest/>`_ for more information.

Installation from source
------------------------
Since |hazel2| is still in a heavy development phase, the only possible way of
installing the code is via compilation of the sources.
Once your environment is activated (or in your base environment if you prefer not to
use a different environment), the sources can be installed by cloning this repository and installing it:

::
    
    git clone https://github.com/aasensio/hazel2
    python setup.py install

Do not forget to often pull from the `repository <https://github.com/aasensio/hazel2>`_ and recompile
the code by typing the following from the location of the sources:

::
    
    git pull
    python setup.py install
    
Improvements are pushed all the time in |hazel2|. In case you do not want
to reinstall it after every pull or you are in fact helping on the development
of the code, you can also install it in develop mode:

::
    
    git pull
    python setup.py develop


..
    Using pip (not preferred until a stable version is published)
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
* ``mpich``
* ``mpi4py``
* ``configobj``
* ``tqdm``
* ``asciitree``

If you want to use the GUI, you also need to install the following dependencies:

* ``matplotlib``
* ``pyqt5``

All of them can be installed in Anaconda with:

::

    conda install numpy h5py scipy astropy mpi4py configobj tqdm pyqt maplotlib

For developers
--------------
Remember that if you want to be involved in the development of the code, it is perhaps more
useful to install the code using

::
    
    python setup.py develop

so that you can immediately test the changes to the code. This development
version can be uninstalled by typing:

::
    
    python setup.py develop --uninstall

Documentation
-------------
For the compilation of the docs, you will need the following packages:

::
  
    sphinx numpydoc sphinx_rtd_theme nbsphinx

Issues
------

MacOSX Catalina has some issues with the paths to the include files. Defining the following environment variables might help.

::

    export CFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
    export CCFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
    export CXXFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
    export CPPFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"

Supercomputers
--------------

Configuring |hazel2| for running in supercomputers should be easy. We recommend
to locally create an Anaconda installation with the needed packages and follow
the details given by the supercomputer staff.
