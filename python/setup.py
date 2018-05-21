# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(
        Extension(
            "DominoChain",
            sources = ["../src/DominoChain.cpp", "./DominoChain.pyx"],
            include_dirs = ["../include", "../../gsl/build-dir",
                            "../../opencv/build/include",
                            numpy.get_include()],
            library_dirs = ["../../gsl/build-dir/Debug",
                            "../../opencv/build/x64/vc15/lib"],
            libraries = ["gsl", "gslcblas", "opencv_world341", "opencv_world341d"],
            language = "c++",
            extra_compile_args = ["-DHAVE_INLINE", "-DVIDEO"]
        )))
