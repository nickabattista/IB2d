from libcpp.list cimport list

cdef extern from "printlist.h":
    void printlist(list[int] &)

def list_test(list[int] l):
    printlist(l)
