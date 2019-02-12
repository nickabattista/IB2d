## Compiling Cython/C data-write methods

Included with the Python code are Python/Cython wrapped C-language data write methods that, when used, improve on the total program running time (exact speedups vary widely based on the simulation - it could be as low as 2% or as high as 100%). Compiling these C methods requires that Cython be installed, which is available as part of the standard Anaconda Python distribution. The libraries will then need to be compiled to run on your computer, after which the Python code should automatically detect them.

In order to compile the wrapped C file-write methods, run the following line in a terminal/command prompt in the pyIB2d/IBM_Blackbox directory:

`python setup.py build_ext —-inplace`

If it succeeds, there will likely be a build directory, a .so or .dll file, and a file wrtFile.c which will appear in the pyIB2d directory.
Then, the IBM_Driver.py will use the C print out method instead of python print out method.

If you are on Windows and get an error message about vcvarsall.bat, try installing `libpython` into your Anaconda environment.
Note however, that getting code to compile for use in Python has been notoriously difficult on Windows and is something of
a moving target. We suggest seeking guidance at https://python-at-risoe.pages.windenergy.dtu.dk/compiling-on-windows/index.html.

Another possible source of error on all machines: calling `python` from a terminal will default to the system python
unless Anaconda python has been prefixed to the shell path (e.g. in .bash_profile - do not prefix user installed
python environments to the system wide path as this could cause major errors!). Make sure you are running the intended
python when you build!!