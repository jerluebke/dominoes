# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(
        Extension(
            "DominoChain",
            sources = ["../src/DominoChain.cpp"],
            include_dirs = ["../include", "../../gsl/build-dir"],
            library_dirs = ["../../gsl/build-dir/Debug"],
            libraries = ["gsl", "gslcblas"]
        ),
        language = "c++11" 
    ))
