from __future__ import print_function
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system
from subprocess import Popen, PIPE
import os

def get_libgfortran_dir():
    """
    Helper function returning the library directory of libgfortran. Useful
    on OSX where the C compiler oftentimes has no knowledge of the library
    directories of the Fortran compiler. I don't think it can do any harm on
    Linux.
    """
    for ending in [".3.dylib", ".dylib", ".3.so", ".so"]:
        try:
            p = Popen(['gfortran', "-print-file-name=libgfortran" + ending],
                      stdout=PIPE, stderr=PIPE)
            p.stderr.close()
            line = p.stdout.readline().decode().strip()
            p.stdout.close()
            if os.path.exists(line):
                return [os.path.dirname(line)]
        except:
            continue
        return []

# compile the fortran modules without linking
# fortran_mod_comp = 'cd src/hazel ; make clean ; make version=python compiler=gfortran'
fortran_mod_comp = 'cd src/hazel ; make version=python compiler=gfortran'
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         "pyhazel",
                         # source file:
                         ['src/hazel/pyhazel.pyx'],
                         libraries=["gfortran"],
                         library_dirs=get_libgfortran_dir(),
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['src/hazel/vars.o', 'src/hazel/maths.o', 'src/hazel/allen.o', 'src/hazel/svd.o', 
                                        'src/hazel/io_py.o', 'src/hazel/SEE.o', 'src/hazel/rt_coef.o', 'src/hazel/synth.o',
														'src/hazel/hazel_py.o','src/hazel/singleton.o'])]

setup(name = 'hazel',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

# system('cp *.so ../pyGUI')
