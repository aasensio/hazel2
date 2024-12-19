#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Hazel v2.0

Hazel synthesis and inversion code. It can invert both photospheric lines and chromospheric lines.

For more information see: https://github.com/aasensio/hazel2
::
    Main Changes in 0.1
    ---------------------
    * Working version
:copyright:
    A. Asensio Ramos    
:license:
    The MIT License (MIT)
"""
import setuptools
from setuptools import find_packages, setup
from setuptools.extension import Extension
from distutils.ccompiler import CCompiler
from distutils.errors import DistutilsExecError, CompileError
from distutils.unixccompiler import UnixCCompiler

from setuptools.command.build_ext import build_ext
from distutils.dep_util import newer_group
from distutils import log
from Cython.Build import cythonize

import os
import platform
from subprocess import Popen, PIPE, check_output
import sys
import numpy
import glob
import re

# Adding an alternative way compatible with PEP440:
def get_pep440version_info():
    # This version of the function uses the format {tag_name}.dev{num_commits}+{commit_hash} to
    # indicate that the version number is a development version, with {tag_name} indicating the
    # name of the latest tag, {num_commits} representing the number of commits since the tag, and
    # {commit_hash} representing the abbreviated commit hash.
    git_describe_output = check_output(['git', 'describe', '--tags']).decode('utf-8').strip()
    
    # Split the output string
    parts = git_describe_output.split('-')
    if len(parts) == 1:
        # No additional commits since the last tag, so this is a release version
        return parts[0]
    else:
        # There are additional commits since the last tag, creating a development version
        tag_name, num_commits, commit_hash = parts
        commit_hash = commit_hash[1:]  # Remove the 'g' prefix from the commit hash
        return f"{tag_name}.dev{num_commits}+{commit_hash}"

DOCSTRING = __doc__.strip().split("\n")

tmp = open('hazel/__init__.py', 'r').read()
author = re.search('__author__ = "([^"]+)"', tmp).group(1)
version = get_pep440version_info()

def _compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
    compiler_so = self.compiler_so    
    arch = platform.architecture()[0].lower()
    if (ext == ".f" or ext == ".f90"):
        if sys.platform == 'darwin' or sys.platform.startswith('linux'):
            compiler_so = [os.environ['FC']]
            if (ext == ".f90"):
                cc_args = ["-O3", "-march=native", "-fPIC", "-c", "-ffree-form", "-ffree-line-length-none", "-w"]
#                cc_args = ["-O3", "-fPIC", "-c", "-ffree-form", "-ffree-line-length-none", "-fno-automatic", "-ffast-math", "-funroll-loops"]
            if (ext == ".f"):
                cc_args = ["-O3", "-march=native", "-fPIC", "-c", "-fno-automatic", "-ffixed-line-length-none", "-w"]
#                cc_args = ["-O3", "-fPIC", "-c", "-ffixed-line-length-none", "-fno-automatic", "-ffast-math", "-funroll-loops"]
            # Force architecture of shared library.
            if arch == "32bit":
                cc_args.append("-m32")
            elif arch == "64bit":
                cc_args.append("-m64")
            else:
                print("\nPlatform has architecture '%s' which is unknown to "
                      "the setup script. Proceed with caution\n" % arch)
    try:        
        self.spawn(compiler_so + cc_args + [src, '-o', obj] +
                   extra_postargs)
    except DistutilsExecError as msg:
        raise CompileError(msg)
UnixCCompiler._compile = _compile


# Hack to prevent build_ext from trying to append "init" to the export symbols.
class finallist(list):
    def append(self, object):        
        return


class MyExtension(Extension):
    def __init__(self, *args, **kwargs):
        Extension.__init__(self, *args, **kwargs)        
        self.export_symbols = finallist(self.export_symbols)        

def get_libgfortran_dir():
    """
    Helper function returning the library directory of libgfortran. Useful
    on OSX where the C compiler oftentimes has no knowledge of the library
    directories of the Fortran compiler. I don't think it can do any harm on
    Linux.
    """
    for ending in [".5.dylib", ".3.dylib", ".dylib", ".3.so", ".so"]:
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

pathGlobal = "src/"

# Monkey patch the compilers to treat Fortran files like C files.
CCompiler.language_map['.f90'] = "c"
UnixCCompiler.src_extensions.append(".f90")
CCompiler.language_map['.f'] = "c"
UnixCCompiler.src_extensions.append(".f")


# SIR
path = pathGlobal+"sir"
list_files = glob.glob(path+'/*.f*')
list_files.append(path+'/sir_code.pyx')

lib_sir = MyExtension('hazel.codes.sir_code',
                  libraries=["gfortran"],
                  library_dirs=get_libgfortran_dir(),
                  sources=list_files,
                  include_dirs=[numpy.get_include()])


# Hazel
path = pathGlobal+"hazel"
list_files = [path+'/vars.f90', path+'/singleton.f90', path+'/maths.f90', path+'/allen.f90', path+'/svd.f90', 
            path+'/io_py.f90', path+'/SEE.f90', path+'/rt_coef.f90', path+'/synth.f90',
			path+'/hazel_py.f90', path+'/hazel_code.pyx']

lib_hazel = MyExtension('hazel.codes.hazel_code',
                  libraries=["gfortran"],
                  library_dirs=get_libgfortran_dir(),
                  sources=list_files,
                  include_dirs=[numpy.get_include()])

class BuildExtSubclass(build_ext):
    def build_extension(self, ext):
        """
        We monkey-patch the build extension function of the build_ext class to avoid sorting.
        This can probably be done in a more elegant way, but it works. Files in F90 cannot 
        be compiled in arbitrary orders because of the module dependencies so the sorted
        introduced in setuptools cannot be used.
        """
        sources = ext.sources
        if sources is None or not isinstance(sources, (list, tuple)):
            raise DistutilsSetupError(
                  "in 'ext_modules' option (extension '%s'), "
                  "'sources' must be present and must be "
                  "a list of source filenames" % ext.name)
        # This sorting needs to be removed
        # sources = sorted(sources)

        ext_path = self.get_ext_fullpath(ext.name)
        depends = sources + ext.depends
        if not (self.force or newer_group(depends, ext_path, 'newer')):
            log.debug("skipping '%s' extension (up-to-date)", ext.name)
            return
        else:
            log.info("building '%s' extension", ext.name)

        # First, scan the sources for SWIG definition files (.i), run
        # SWIG on 'em to create .c files, and modify the sources list
        # accordingly.
        sources = self.swig_sources(sources, ext)

        # Two possible sources for extra compiler arguments:
        #   - 'extra_compile_args' in Extension object
        #   - CFLAGS environment variable (not particularly
        #     elegant, but people seem to expect it and I
        #     guess it's useful)
        # The environment variable should take precedence, and
        # any sensible compiler will give precedence to later
        # command line args.  Hence we combine them in order:
        extra_args = ext.extra_compile_args or []

        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))

        objects = self.compiler.compile(sources,
                                         output_dir=self.build_temp,
                                         macros=macros,
                                         include_dirs=ext.include_dirs,
                                         debug=self.debug,
                                         extra_postargs=extra_args,
                                         depends=ext.depends)

        # XXX outdated variable, kept here in case third-part code
        # needs it.
        self._built_objects = objects[:]

        # Now link the object files together into a "shared object" --
        # of course, first we have to figure out all the other things
        # that go into the mix.
        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []

        # Detect target language, if not provided
        language = ext.language or self.compiler.detect_language(sources)

        self.compiler.link_shared_object(
            objects, ext_path,
            libraries=self.get_libraries(ext),
            library_dirs=ext.library_dirs,
            runtime_library_dirs=ext.runtime_library_dirs,
            extra_postargs=extra_args,
            export_symbols=self.get_export_symbols(ext),
            debug=self.debug,
            build_temp=self.build_temp,
            target_lang=language)


    #     arch = platform.architecture()[0].lower()
    #     if (ext == ".f" or ext == ".f90"):
    #         if sys.platform == 'darwin' or sys.platform.startswith('linux'):
    #             compiler_so = [os.environ['FC']]
    #             if (ext == ".f90"):
    #                 cc_args = ["-O3", "-march=native", "-fPIC", "-c", "-ffree-form", "-ffree-line-length-none", "-w"]
    # #                cc_args = ["-O3", "-fPIC", "-c", "-ffree-form", "-ffree-line-length-none", "-fno-automatic", "-ffast-math", "-funroll-loops"]
    #             if (ext == ".f"):
    #                 cc_args = ["-O3", "-march=native", "-fPIC", "-c", "-fno-automatic", "-ffixed-line-length-none", "-w"]
    # #                cc_args = ["-O3", "-fPIC", "-c", "-ffixed-line-length-none", "-fno-automatic", "-ffast-math", "-funroll-loops"]
    #             # Force architecture of shared library.
    #             if arch == "32bit":
    #                 cc_args.append("-m32")
    #             elif arch == "64bit":
    #                 cc_args.append("-m64")
    #             else:
    #                 print("\nPlatform has architecture '%s' which is unknown to "
    #                     "the setup script. Proceed with caution\n" % arch)
    #     breakpoint()

### https://github.com/manodeep/Makefile-C-Python/blob/master/setup.py

setup_config = dict(
    name='hazel',
    version=version,
    description=DOCSTRING[0],
    long_description="\n".join(DOCSTRING[2:]),
    author=author,
    author_email='aasensio@iac.es',
    url='https://github.com/aasensio/hazel2',
    license='GNU General Public License, version 3 (GPLv3)',
    platforms='OS Independent',
    install_requires=['numpy','scipy','configobj','h5py','astropy','tqdm','cython'],
    ext_modules=cythonize([lib_sir, lib_hazel]),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: Implementation :: CPython",
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=['hazel', 'radiative transfer'],
    packages=find_packages(),
    zip_safe=False,
    include_package_data=True,
    scripts=['gui/hazelgui'],
    cmdclass={'build_ext': BuildExtSubclass}
)

if __name__ == "__main__":    
    setup(**setup_config)

#    Attempt to remove the *.c files
    for filename in glob.glob("**/*.c", recursive=True):
        try:
            os.remove(filename)
        except:
            pass

    for filename in glob.glob("*.mod", recursive=True):
        try:
            os.remove(filename)
        except:
            pass
