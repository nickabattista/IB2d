#!/usr/bin/env python

'''To compile: python setup.py build_ext --inplace'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("write",
                             sources=["write.pyx", "c_write.c"],
                             include_dirs=[numpy.get_include()])],
)
