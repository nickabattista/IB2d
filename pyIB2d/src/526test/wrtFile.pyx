"""
wrtFile.pyx

simple cython test of accessing a numpy array's data

the C function: c_multiply multiplies all the values in a 2-d array by a scalar, in place.

"""

import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_savevtk_points_write(int N, double *x, char* filename, char* vectorName)

@cython.boundscheck(False)
@cython.wraparound(False)
def write(int n,np.ndarray[double, ndim=2, mode="c"] input not None,str filename, str vecname):
    filename_by = filename.encode('UTF-8')
    vecname_by = vecname.encode('UTF-8')
    cdef char* c_filename = filename_by
    cdef char* c_vecname = vecname_by
    c_savevtk_points_write(n,&input[0,0], c_filename, c_vecname)

    return None