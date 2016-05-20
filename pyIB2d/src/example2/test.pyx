cimport numpy as np
import numpy as np

np.import_array()

cdef extern from "printlist.h":
    void printlist(list[int] &)

def list_test(list[int] l):
    printlist(l)
