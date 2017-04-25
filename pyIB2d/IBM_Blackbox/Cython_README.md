## Compiling Cython/C data-write methods

Included with the Python code are Python/Cython wrapped C-language data write methods that, when used, improve on the total program running time (exact speedups vary widely based on the simulation - it could be as low as 2% or as high as 100%). Compiling these C methods requires that Cython be installed, which is available as part of the standard Anaconda Python distribution. The libraries will then need to be compiled to run on your computer, after which the Python code should automatically detect them.

In order to compile the wrapped C file-write methods, run the following line in a terminal/command prompt in the pyIB2d directory:

python setup.py build_ext —-inplace

If it succeeds, there will likely be a build directory, a .so or .dll file, and a file wrtFile.c which will appear in the pyIB2d directory.
Then, the IBM_Driver.py will use the C print out method instead of python print out method.
