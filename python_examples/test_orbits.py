#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
import time
import copy as C

from mpl_toolkits.mplot3d import Axes3D

# Tracking
import Trpy as Tr

# detectors
import detectors as Det

# limiter drawing
import get_limiter as gl

dtr = np.pi/180.

det_colors = ['r','g','b', 'y', 'm', 'c', 'o', 'k']

#%%  set initial valuesfor tracker
# 
Tr.tracker.particle_charge = 1.
# Tr.tracker.particle_mass_amu = 1.007347  
Tr.tracker.particle_mass_amu = 1.00727647   # this is not really accurate
Tr.tracker.particle_energy_mev = 3.
Tr.tracker.track_length = 6.
Tr.tracker.step_size = 0.01
Tr.control_mod.time_reversed = True

# Select which tracker to use
#Tr.tracker.selected_tracker = Tr.tracker.bulirsch_stoer_t
Tr.tracker.selected_tracker = Tr.tracker.boris_t

# file names should be set using these functions
Tr.tracker.set_gfile_name('029880.00234.dat')
Tr.tracker.set_efit_directory('../example_data/')

#%% initialize tracker: allocate space for trajectory data in BT.tracker.trajectory
# 
Tr.tracker.init_tracker()
#%% load flux i.e the magnetic field and its interpolations
# 
Tr.tracker.load_flux()


#%% setup the limiter

Tr.limiter_control_mod.set_limiter_file_name('MASTLIMIT00.dat')
Tr.limiter_control_mod.set_limiter_directory('../example_data/')
Tr.limiter_control_mod.limiter_init()
# setup for limiter drawing
limiter = gl.limiter('limiter_drawing.data')
Tr.control_mod.check_limiter = True
Tr.control_mod.print_hit = False

#%% setup a detector

det1 = Det.detector(1, 'Detector 1', 
                    position = np.array([1.6, 0., 0.*dtr]), 
                    direction = np.array([45*dtr, 0.*dtr]), 
                    rotation = 0.*dtr,
                    tracker = Tr,
                    bundle_fname = 'Detector1.npz')
# set b-field 0 for test
# det1.B_det = np.array([0.,0.,0.])
det1.init_trajectories()
det1.calc_trajectories()

det2 = Det.detector(2, 'Detector 2', 
                    position = np.array([1.6, 0., 0.*dtr]),
                    pos_local = np.array([0.05, 0., 0.05]),
                    direction = np.array([45*dtr, 0.*dtr]), 
                    rotation = 0.*dtr,
                    tracker = Tr,
                    bundle_fname = 'Detector2.npz')
# set b-field 0 for test
# det1.B_det = np.array([0.,0.,0.])
det2.init_trajectories()
det2.calc_trajectories()

det_l = [det1, det2]

#%% Plotting
# close all figures using close('all') 
# make 3d plot
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

# plot all  in det_l
for dd in det_l:
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
for dd in det_l:
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
for dd in det_l:
    for b in dd.bundle:
        x = b[0]
        y = b[1]
        B.pl.plot(x,y, color = det_colors[dd.number])
B.pl.ylim((-2.,2.))
B.pl.xlim((-2.,2.))

B.pl.show()


