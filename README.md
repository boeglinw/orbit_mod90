# orbit_mod90
Modular Lorentz Orbit Code with Python Wrappers. This is a new version of orbit3. For tracking it uses the Boris algorithm. 

## installation:
Run the Makefile in the root directory. This will compile the libraries and create the python wraper using f2py which needs 
to be installes. 

##Examples
In the directory *example_data* are data examples and python scripts to explore and test the modules. The python modules are located in the *python* directory.
This includes the shared library BTpy.so (after the make) and a few example python scripts. In various directories are also fortan90 example codes to test the
different libraries. All has been compiled with gfortran (GNU Fortran (GCC) 8.2.0) and f2py (Version: 2, numpy Version: 1.18.1)

