Wrapped C print methods are provided that will approximately give a 6% improvement on the total program running time, thanks to the fast speed of processing character arrays. Compiling these C methods requires Cython, which is available as part of the standard Anaconda Python distribution.

In order to use the wrapped C print-out methods, run the following line in a terminal/command prompt in the pyIB2d directory:

python setup.py build_ext —-inplace

If it succeeds, there will be a build directory, a .so or .dll file, and a file wrtFile.c which will appear in the pyIB2d directory.
Then, the IBM_Driver.py will automatically use the C print out method instead of python print out method.
