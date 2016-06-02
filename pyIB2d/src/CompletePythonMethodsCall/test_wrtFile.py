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

colorMap = 'RED'

dx = 3.0
dy = 4.0

write.savevtk_points_write(Arow, A,'savevtk_points_write.txt', vecname)
write.savevtk_points_connects_write(Arow, Brow, A,'savevtk_points_connects_write.txt',vecname,B)
write.savevtk_vector(Arow,Acol,Arow,Acol, A, A,'savevtk_vector.txt',vecname,dx,dy)
write.savevtk_scalar(Arow,Acol,A,'savevtk_scalar.txt', colorMap,dx,dy)