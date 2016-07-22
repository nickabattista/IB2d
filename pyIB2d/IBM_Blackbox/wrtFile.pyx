"""
wrtFile.pyx

simple cython test of accessing a numpy array's data

the C function: c_multiply multiplies all the values in a 2-d array by a scalar, in place.

Created by Ao Zeng on 5/17/16.
"""
import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_savevtk_points_write(int N, double *x, char* filename, char* vectorName)
cdef extern void c_savevtk_points_connects_write(int N, int Nc, double* x,char* filename, char* vectorname, double* connectsMat)
cdef extern void c_savevtk_vector(int Xrow, int Xcol, int Yrow, int Ycol, double* X, double* Y,char* filename,char* vectorName,double dx,double dy)
cdef extern void c_savevtk_scalar(int arrayRow,int arrayCol,double* array, char* filename, char* colorMap,double dx, double dy)

@cython.boundscheck(False)
@cython.wraparound(False)
def savevtk_points_write(int n,np.ndarray[double, ndim=2, mode="c"] input not None,str filename, str vecname):
    filename_by = filename.encode('UTF-8')
    vecname_by = vecname.encode('UTF-8')
    cdef char* c_filename = filename_by
    cdef char* c_vecname = vecname_by
    c_savevtk_points_write(n,&input[0,0], c_filename, c_vecname)

    return None

def savevtk_points_connects_write(int n, int nc, np.ndarray[double, ndim=2, mode="c"] x,str filename, str vecname, np.ndarray[double, ndim=2, mode="c"] connectsMat):
    filename_by = filename.encode('UTF-8')
    vecname_by = vecname.encode('UTF-8')
    cdef char* c_filename = filename_by
    cdef char* c_vecname = vecname_by
    c_savevtk_points_connects_write(n,nc,&x[0,0],c_filename,c_vecname,&connectsMat[0,0])
    
    return None

def savevtk_vector(int xRow, int xCol, int yRow, int yCol, np.ndarray[double, ndim=2, mode="c"] x, np.ndarray[double, ndim=2, mode="c"] y,str filename,str vecname,double dx,double dy):
    filename_by = filename.encode('UTF-8')
    vecname_by = vecname.encode('UTF-8')
    cdef char* c_filename = filename_by
    cdef char* c_vecname = vecname_by
    c_savevtk_vector(xRow, xCol, yRow, yCol, &x[0,0], &y[0,0],c_filename,c_vecname,dx,dy)
    
    return None

def savevtk_scalar(int arrayRow,int arrayCol,np.ndarray[double, ndim=2, mode="c"] array, str filename, str colorMap,double dx, double dy):
    filename_by = filename.encode('UTF-8')
    colorMap_by = colorMap.encode('UTF-8')
    cdef char* c_filename = filename_by
    cdef char* c_colorMap = colorMap_by
    c_savevtk_scalar(arrayRow,arrayCol,&array[0,0],c_filename, c_colorMap,dx, dy)

    return None







