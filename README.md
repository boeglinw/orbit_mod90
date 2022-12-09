# orbit_mod90
Modular Lorentz Orbit Code with Python Wrappers. This is a new version of orbit3. For tracking it uses the Boris algorithm or a Bulirsch-Stoer and Runge-Kutta integrator. This can be selected the control_module.

## Installation:
The main Makefile will create all the modules and libraries.
Run the Makefile in the root directory as follows:

- make  

This will compile all the required libraries but not the python modules.

To create the python modules you need f2py.
The python modules are then created by doing (again in the root directory):

- make modules



## Examples:
In the directory *example_data* are data examples to explore and test the modules. The python modules are located in the *python_modules* directory.
The following two main python examples are in the directory *python_examples*.

- example_Trpy.py

This program calculates trajectories in the equilibrium magnetic field starting at a certain initial position and velocity.

- calc_orbits.py

This is a more complete program that calculates orbit bundles for a set of cylindrical detectors. The necessary parameters are controlled by the *calc\_orbit\_control.data* file.

All has been compiled with *gfortran (GNU Fortran (GCC) 8.2.0)* and *f2py (Version: 2, numpy Version: 1.19.2)*
Most of the example scripts require *LabTools3* which can be installed using *pip (pip install LabTools3)* or from GitHub: *https://github.com/boeglinw/LabTools3.git*

It is recommended to calculate orbits in a working directory instead of the orbit_mod90 directory. You can copy example parameter files and the calc_orbits.py file there. Please also read the top of the comments in calc_orbits.py about the setting of the PYTHONPATH to be able to run the code.