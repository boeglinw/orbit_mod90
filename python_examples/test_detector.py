#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 07:52:36 2020

@author: boeglinw
"""
import calc as C
import numpy as np
import matplotlib.pyplot as pl

# Tracking
import Trpy as Tr
# detectors
import detectors as Det

import acceptance_array as AA

import copy as CO

from numba import jit

# timing purposes
import timeit

import fibonacci as FI
dtr = np.pi/180.

# initialize
fg = FI.fibonacci_grid()

#%% kinematics
# kinetic energy

K = 3*C.mev

# velocity
v = np.sqrt(2.*K/C.mp)

mp_me = C.mp/C.me

"""+
thv = 0.*dtr
phv = 0.*dtr

vs = np.array([np.sin(thv)*np.cos(phv),
               np.sin(thv)*np.sin(phv),
               np.cos(thv)])*v

rs = np.array([-0.0,0.,0.])
"""
Nt = 50   # number of steps in straight section

Rd = 0.0019
Rc = 0.0019
Rcyl = 0.0019
d = 0.0357

step = d/Nt


# conical detector
h_c = AA.con_detector(Rd, Rc, d)
h_c.calc_acceptance()

# largest angle
th_max = np.arctan((Rc + Rd)/d)*1.
# th_max = np.pi/2.5

# set  B-field
Bmag = 1.0
#Bmag = 0.

theta_b = 90*dtr
phi_b = 0.*dtr

# magnetic field direction in detector system
Bv = np.array([np.cos(phi_b)*np.sin(theta_b), np.sin(phi_b)*np.sin(theta_b), np.cos(theta_b)])*Bmag

#%%
# set initial valuesfor tracker
Tr.tracker.particle_charge = 1.
# Tr.tracker.particle_mass_amu = 1.007347
Tr.tracker.particle_mass_amu = 1.00   # this is not really accurate
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

#%% setup detector
# sunflower r theta distribution for initial position
np_pos= 50 # set the number of points

# sunflower for directions

th_scale = 3
np_dir = int(100*th_scale**2)

R_p, Z_p, Phi_p = 1.639, 0., 80.3656*dtr 

xd, yd, zd = 0., -28.11e-3, 24.414e-3 

theta_d, phi_d = 40.*dtr, -45.*dtr 

det_l = Det.detector(0,
            'Channel_0',
            head_position = np.array([R_p, Z_p, Phi_p]),
            det_pos_local = np.array([xd, yd, zd]),
            det_direction = np.array([theta_d, phi_d]),
            arm_rotation = 0.,
            tracker = Tr,
            bundle_fname = 'Channel_0_test.npz',
            color = 'red',
            R_det = Rd, R_coll = Rc, R_cyl = Rc, D = d, 
            zero_at_coll = False)



#%% calculations
# init
det_l.init_trajectories(N_pos = 20, N_dir = 100., N_s = Nt, dN_s = 5)
det_l.calc_trajectories()

"""
#%% analze results

# compare to analytical calc.
Ratio = bd.Acc/h_c.accept

print(f'Ratio = {bd.Acc/h_c.accept}')



#%% plot final positions
pl.figure()
pl.plot(r_f[0], r_f[1], '.', color = 'b')
pl.gca().set_aspect('equal')

#%%
cth = v_f[2]/v
"""


