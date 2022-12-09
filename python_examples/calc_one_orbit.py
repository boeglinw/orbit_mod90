#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

An esxample if using the tracker to calculate a single trjectory

One can select initial position and velocity and either select time reversed or

normal calculation (time reversed changes the direction of the magntic field)

@author: boeglinw
"""

import numpy as np
import LT.box as B
import time
import sys

import copy as C

from mpl_toolkits.mplot3d import Axes3D

#%% setup the PYTHON path

def add_sys_path(new_path):
    """ AddSysPath(new_path): adds a directory to Python's sys.path

    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os

    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1

    # Standardize the path. Windows is case-insensitive, so lowercase
    # for definiteness.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower(  )

    # Check against all currently available paths
    for x in sys.path:
        x = os.path.abspath(x)
        if sys.platform == 'win32':
            x = x.lower(  )
        if new_path in (x, x + os.sep):
            return 0
    sys.path.append(new_path)
    return 1

#----------------------------------------------------------------------------------------------------------
#%% location of the python modules for orbit_mod90 : Adjus this to your system
orbit_mod90_python = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/python_modules'
#----------------------------------------------------------------------------------------------------------


if (add_sys_path(orbit_mod90_python)) < 0 :
    print(f'Cannot add {orbit_mod90_python} to PYTHONPATH as it does not exist !')
    sys.exit(-1)

#%% import Trpy
import Trpy as Tr
import coordinate_systems as Cs

#%%  some helper functions
# for simulation purposes only
twopi = 2.*np.pi
def pol_angle(x,y):
    # get phi for a vector
    angle = (np.arctan2(y,x) + twopi)%twopi
    return angle

# get uniformly distributed pints inside a circle of radius r
def get_random_point(r):
    outside = True
    while outside:
       x = 2.*(np.random.uniform() - 0.5)*r
       y = 2.*(np.random.uniform() - 0.5)*r
       outside = (x**2 + y**2 > r)
    return x,y


#%% setup kinematics and tracking parameters


Tr.tracker.particle_charge_ec = 1.
# Tr.tracker.particle_mass_amu = 1.007347
Tr.tracker.particle_mass_amu = 1.00   # this is not really accurate
Tr.tracker.particle_energy_mev = .0008
Tr.tracker.track_length = 100.
Tr.tracker.step_size = 0.002
Tr.control_mod.time_reversed = False  # ths is not working correctly needs to be checked

# Select which tracker to use
#Tr.tracker.selected_tracker = Tr.tracker.bulirsch_stoer_t
Tr.tracker.selected_tracker = Tr.tracker.boris_t

#%% get magntic field

# file names should be set using these functions
Tr.tracker.set_gfile_name('029880.00234.dat')
Tr.tracker.set_efit_directory('../example_data/')

#%%
# initialize tracker: allocate space for trajectory data in BT.tracker.trajectory
Tr.tracker.init_tracker()

#%%
# load flux i.e the magneit field and its interpolations
Tr.tracker.load_flux()

#%% setup the limiter

Tr.limiter_control_mod.set_limiter_file_name('MASTLIMIT00.dat')
Tr.limiter_control_mod.set_limiter_directory('../example_data/')
Tr.limiter_control_mod.limiter_init()
Tr.control_mod.check_limiter = True
Tr.control_mod.print_hit = True

#%%
# set initial track position

# enpoint of a 3 MeV trajectory starting at [1., 0., 0.] in direction [0.,0.,1]
# r0 = np.array([1.9472439, 0.4492197, 0.6346036])
# this would be the initial velocity
# vi = -np.array([20381802.877327, 12203044.718455, -3818356.074167])

# final velocity

# some random position inside the plasma
r0 = np.array([1.2, 0., 0.])


# initial position in R-Z
R0 = np.sqrt(r0[0]**2 + r0[1]**2)
Z0 = r0[2]




#%% magnetic field at initial position

# for starting partivles in a direction close to B
b_pol_r, b_pol_z, b_phi, b_total, psi_rel = Tr.em_fields_mod.bfield(R0,Z0)

r0_mag = Cs.L.norm(r0)

er = np.array([r0[0]/r0_mag, r0[1]/r0_mag, 0.])
ephi = np.array([-er[1], er[0], 0.] )
ez = np.array([0., 0., 1.])

Bf0 = b_pol_r*er + b_pol_z*ez + b_phi*ephi

uBf0 = Cs.get_unit_vector(Bf0)

# perpendicular perturbation
dBf0 = np.cross(r0, Bf0)
udBf0 = Cs.get_unit_vector(dBf0)

# size pf perpendicular compilnent compared to 1 for the paralell
eps = 1.1031


#%% initial velocity unit vector

# various optyionms for the initial velocity in the normal coordinate system

# vi = np.array([0,0,1.])


#vi = -np.array([20381802.877327, 12203044.718455, -3818356.074167])


vi = Cs.get_unit_vector(uBf0 + eps*udBf0)


uv0 = Cs.get_unit_vector( vi )
# set initial velocity
v0 = Tr.tracker.vmag*uv0

#%%
t_start = time.time()  # for timing

# calculate track result is in Tt.tracker.trajectory array
nc = Tr.tracker.get_trajectory(r0, v0 )
print ('trajectory contains ', nc, ' steps' )

# need to copy the result values to store them
x = C.copy(Tr.tracker.trajectory[:nc-1,0])
y = C.copy(Tr.tracker.trajectory[:nc-1,1])
z = C.copy(Tr.tracker.trajectory[:nc-1,2])

vx = C.copy(Tr.tracker.trajectory[:nc-1,3])
vy = C.copy(Tr.tracker.trajectory[:nc-1,4])
vz = C.copy(Tr.tracker.trajectory[:nc-1,5])

r = np.sqrt(x**2 + y**2)

bf = np.array([Tr.em_fields_mod.bfield(rr,zz) for rr,zz in zip(r,z)])


#%% Plotting
# close all figures using close('all')
# make 3d plot
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

B.pl.plot(x, y, z, color = 'r')
# B.pl.plot(xpz,ypz,zpz, color = 'g')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


#%% compare 2d plot
fig2 = B.pl.figure()
ax2 = fig2.add_subplot(111)
B.pl.plot(r,z, color = 'r')
ax2.set_aspect('equal')
ax2.set_xlabel('R')
ax2.set_ylabel('Z')



