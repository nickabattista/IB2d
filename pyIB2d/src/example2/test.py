#!/usr/bin/env python

"""
simple test of the multiply.pyx and c_multiply.c test code
"""


import numpy as np

import write

a = np.array([[1.0, 2.0], [3.0, 4.0]])

print (a)

write.write(a)

print (a)