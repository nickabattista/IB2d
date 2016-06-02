#!/usr/bin/env python3

'''Python file to test c_wrtFile.c and wrtFile.pyx'''

import numpy as np
import write

A = np.array([[1.,2.,3.],[4.,5.,6.]])
Arow = A.shape[0]
Acol = A.shape[1]

B = np.random.rand(2,3)
Brow = B.shape[0]
Bcol = B.shape[1]

filename = 'itworks.txt'

vecname = 'foo'

write.write(Arow, A,filename, vecname)
