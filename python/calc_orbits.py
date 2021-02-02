#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

Main driving script to calculate orbits or orbit bundles

controlled by input file

@author: boeglinw
"""

import numpy as np
import LT.box as B
from LT.parameterfile import pfile

import sys


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


#%% location of the python modules for orbit_mod90 
orbit_mod90_python = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/python'

if (add_sys_path(orbit_mod90_python)) < 0 :
    print(f'Cannot add {orbit_mod90_python} to PYTHONPATH as it does not exist !')
    sys.exit(-1)


# Tracking
import Trpy as Tr

# detectors
import detectors as Det

# limiter drawing
import get_limiter as gl

dtr = np.pi/180.

det_colors = ['r','g','b', 'y', 'm', 'c', 'o', 'k']

trackers = {'Boris':Tr.tracker.boris_t, 'Bulirsch_stoer':Tr.tracker.bulirsch_stoer_t}


#%%  control file input

# load control file, this is a pure parameter file
cd = pfile('example_control.data')

# this is a data file describing the detector head
dh = B.get_file(cd['detector_head'])


#%%  set initial valuesfor tracker
# 
Tr.tracker.particle_charge = cd.get_value('particle_charge')
# Tr.tracker.particle_mass_amu = 1.007347  
Tr.tracker.particle_mass_amu = cd.get_value('particle_mass_amu')
Tr.tracker.particle_energy_mev = cd.get_value('particle_energy_mev')

Tr.tracker.track_length = cd.get_value('track_length')
Tr.tracker.step_size = cd.get_value('step_size')
Tr.control_mod.time_reversed = cd.get_value('time_reversed')

# Select which tracker to use
Tr.tracker.selected_tracker = trackers[cd.get_value('selected_tracker')]

# file names should be set using these functions
Tr.tracker.set_gfile_name(cd.get_value('gfile_name'))
Tr.tracker.set_efit_directory(cd.get_value('efit_directory'))


# Control trajectory bundle calculation
# number of fibonacci positions in detector area
N_pos_det = cd.get_value('number_of_fib_points')
# number of fibonacci points for direction
N_dir_det = cd.get_value('number_of_fib_directions')


#%% initialize tracker: allocate space for trajectory data in BT.tracker.trajectory
# 
Tr.tracker.init_tracker()
#%% load flux i.e the magnetic field and its interpolations
# 
Tr.tracker.load_flux()


#%% setup the limiter

Tr.limiter_control_mod.set_limiter_file_name(cd.get_value('limiter_file_name'))
Tr.limiter_control_mod.set_limiter_directory(cd.get_value('limiter_directory'))
Tr.limiter_control_mod.limiter_init()
Tr.control_mod.check_limiter = cd.get_value('check_limiter')
Tr.control_mod.print_hit = cd.get_value('print_hit')
# setup for limiter drawing
limiter = gl.limiter('limiter_drawing.data')
# Tr.control_mod.print_polygon = False

#%% setup a detectors
R_p = cd.get_value('detector_head_R')
Z_p = cd.get_value('detector_head_Z')
Phi_p = cd.get_value('detector_head_Phi')*dtr

arm_rotation = cd.get_value('arm_rotation')*dtr

detector_head = []

for i,n  in enumerate(dh['Detector_number']):
    det_l = Det.detector(n, 
                        dh['Detector_name'][i], 
                        position = np.array([R_p, Z_p, Phi_p]), 
                        pos_local = np.array([dh['xd'][i], dh['yd'][i], dh['zd'][i]]),
                        direction = np.array([dh['theta_d'][i]*dtr, dh['phi_d'][i]*dtr]), 
                        rotation = arm_rotation,
                        tracker = Tr,
                        bundle_fname = f'det_{n}_'+dh['Detector_name'][i]+'.npz')
    detector_head.append(det_l)
# initialize
for det_l in detector_head:
    det_l.init_trajectories(N_pos = N_pos_det, N_dir = N_dir_det)
    
# calculate trajectories
for det_l in detector_head:
    det_l.calc_trajectories()


#%% Plotting
# close all figures using close('all') 
# make 3d plot
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

# plot all  in det_l
for dd in detector_head:
    for b in dd.bundle:
        x = b[0]
        y = b[1]
        z = b[2]
        B.pl.plot(x,y,z, color = det_colors[dd.number])
# B.pl.plot(xpz,ypz,zpz, color = 'g')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


#%% RZ plot
fig_rz = B.pl.figure()
limiter.draw_side_all()

# plot all  in det_l
for dd in detector_head:
    for b in dd.bundle:
        r = b[3]
        z = b[2]
        B.pl.plot(r,z, color = det_colors[dd.number])
B.pl.ylim((-2.,2.))
B.pl.xlim((0.,2.))


#%% midplane plot
fig_mid = B.pl.figure()
limiter.draw_top_all()
# plot all  in det_l
for dd in detector_head:
    for b in dd.bundle:
        x = b[0]
        y = b[1]
        B.pl.plot(x,y, color = det_colors[dd.number])
B.pl.ylim((-2.,2.))
B.pl.xlim((-2.,2.))

B.pl.show()


