.. include:: ../hazel_name

Installation
============

Using pip
---------

|hazel2| is distributed with precompiled binaries for Linux (x86_64 and aarch64) and macOS
(Apple Silicon and Intel), so no compiler is needed. It requires Python 3.11 or newer. We recommend
installing it in a virtual environment:

::

    python -m venv hazel_env
    source hazel_env/bin/activate
    pip install hazel2

Note that the package is called ``hazel2`` on PyPI, but it is imported as ``hazel``. Windows
users can install it inside `WSL <https://learn.microsoft.com/windows/wsl/>`_.

The following optional features can be added between brackets, for instance ``pip install "hazel2[nn,mpi]"``:

* ``nn``: PyTorch, needed for the NLTE neural option for synthesizing the Ca II 854.2 nm line.
* ``mpi``: ``mpi4py``, for inverting large maps in parallel. It needs an MPI library in your system.
* ``zarr``: support for ``zarr`` input and output files.
* ``gui``: ``matplotlib`` and ``PyQt5``, needed for the ``hazelgui`` graphical interface.

The first time the NLTE option is used, the weights of the neural network (76 MB) are downloaded
and stored in ``~/.cache/hazel``. If the nodes of your supercomputer have no access
to the network, download them from a login node:

::

    python -c "import hazel; hazel.download_nlte_model()"

The file can also be placed in any other directory, indicated with the ``HAZEL_DATA_DIR`` environment variable.

Installation from source
------------------------

Since improvements are pushed all the time, you might want to use the development version. This
needs a Fortran and a C compiler (``gfortran`` and ``gcc``, for instance from
``conda install -c conda-forge compilers``). The remaining build tools are downloaded by ``pip``:

::

    git clone https://github.com/aasensio/hazel2
    cd hazel2
    pip install .

Binaries are built to run in any computer. To tune the code for your own CPU, which is a
few percent faster, use:

::

    pip install . -Csetup-args=-Dnative=true

Do not forget to often pull from the `repository <https://github.com/aasensio/hazel2>`_ and reinstall
the code by typing the following from the location of the sources:

::

    git pull
    pip install .

If you want a verbose installation, you can also add the ``-v`` flag to the previous command.

Container
---------

A simple way of installing |hazel2|, at least in Linux systems or supercomputers, is to use 
this Apptainer `container <https://cloud.iac.es/index.php/s/BqG2jGRsqHtJMW9>`_.
This container has all the dependencies installed and it is ready to be used. You can download it (we try to
keep this version as updated as possible) and run it with Apptainer as:

::    

    apptainer run hazel.sif python invert.py

In principle, you can build your own container by going to the `apptainer` directory in the
distribution and running:

::

    apptainer build hazel.sif hazel.def


For developers
--------------
If you want to be involved in the development of the code, install it in editable mode,
so that changes are used immediately and the Fortran code is recompiled when importing ``hazel``:

::

    pip install meson-python "meson<1.11" ninja cython numpy
    pip install --no-build-isolation -e .

This version can be uninstalled by typing:

::

    pip uninstall hazel2

Documentation
-------------
For the compilation of the docs, you will need the following packages:

::
  
    conda install sphinx numpydoc sphinx_rtd_theme nbsphinx pandoc

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
