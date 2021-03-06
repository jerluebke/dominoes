# -*- coding: utf-8 -*-

import sys
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

MAKROS = ["-DHAVE_INLINE", "-DVIDEO"]

if "--print_extra" in sys.argv:
    print("Debug mode ...")
    MAKROS.extend(["-DPRINT_EXTRA"])
    sys.argv.remove("--print_extra")
elif "--rmod" in sys.argv:
    print("using modified R value ...")
    MAKROS.extend(["-DPRINT_EXTRA", "-DR_ORIG=1"])
    sys.argv.remove("--rmod")

setup(ext_modules = cythonize(
        Extension(
            "DominoChain",
            sources = ["../src/DominoChain.cpp", "./DominoChain.pyx"],
            include_dirs = ["../include", "../../gsl/build-dir",
                            "../../opencv/build/include",
                            numpy.get_include()],
            library_dirs = ["../../gsl/build-dir/Debug",
                            "../../opencv/build/x64/vc15/lib"],
            libraries = ["gsl", "gslcblas", "opencv_world341",
                         "opencv_world341d"],
            language = "c++",
            extra_compile_args = MAKROS
        )))
