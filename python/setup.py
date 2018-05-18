# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(
        Extension(
            "DominoChain",
            sources = ["../src/DominoChain.cpp", "./DominoChain.pyx"],
            include_dirs = ["../include", "../../gsl/build-dir",
                            numpy.get_include()],
            library_dirs = ["../../gsl/build-dir/Debug"],
            libraries = ["gsl", "gslcblas"],
            language = "c++"
        )))
