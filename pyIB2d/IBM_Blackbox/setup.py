#!/usr/bin/env python3

'''To compile: python setup.py build_ext --inplace'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("write",
                             sources=["wrtFile.pyx", "c_write.c"],
                             include_dirs=[numpy.get_include()])],
    python_requires='>=3.5',
    )