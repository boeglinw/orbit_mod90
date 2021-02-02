#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
import time
import Trpy as Tr
import copy as C

from mpl_toolkits.mplot3d import Axes3D



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

#%% setup    
track_dir = '../example_data/'
track_name = 'track_21111.data'
# for comparison
td = B.get_file(track_dir + track_name)
rt = td['r']
phit = td['phi']
zt = td['z']
xt = rt*np.cos(phit)
yt = rt*np.sin(phit)

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
Tr.tracker.selected_tracker = Tr.tracker.bulirsch_stoer_t
#Tr.tracker.selected_tracker = Tr.tracker.boris_t

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
r0 = np.array([0.29375441571112459,
               1.6341828183759661,    
               3.8310999999999998E-002])

# initial velocity unit vectopr
uv0 = np.array([0.47137247366158763,
                -0.61971211944698446,      
                0.62750687652382131])

#%% sperical angles of initial velocity
# this is for the generation of a bundle, calculate a trasformation matrix that moves
# the original z-direction into the direction of the initial velocity
phi = pol_angle(uv0[0],uv0[1])
theta = np.arccos(uv0[2])
# rotation matrices
Rz = np.array([ [ np.cos(phi), -np.sin(phi), 0.],
                 [ np.sin(phi), np.cos(phi),  0.],
                 [ 0.,          0.,           1.] ])

Ry = np.array([ [ np.cos(theta),  0., np.sin(theta)],
                 [ 0.,           1.,          0.],
                 [ -np.sin(theta), 0., np.cos(theta)] ])

# rotates the z-axis into the direction of uv0
R_tot = np.matmul(Rz,Ry)

# set initial velocity
v0 = Tr.tracker.vmag*uv0 



#%%
t_start = time.time()  # for timing
# numer of trajectories
n_traj = 10
# store results in bundle
bundle = []
for i in range(n_traj):
    # select a small random offset
    xr,yr = get_random_point(.05)
    dv = np.matmul(R_tot, np.array([xr, yr, 0.]) ) * Tr.tracker.vmag      
    vi = v0 + dv
    # calculate track result is in Tt.tracker.trajectory array
    nc = Tr.tracker.get_trajectory(r0, vi )
    print ('trajectory contains ', nc, ' steps' )
    # need to copy the result values to store them 
    x = C.copy(Tr.tracker.trajectory[:nc-1,0]) 
    y = C.copy(Tr.tracker.trajectory[:nc-1,1])
    z = C.copy(Tr.tracker.trajectory[:nc-1,2])
    r = np.sqrt(x**2 + y**2)
        
    bundle.append((x, y, z, r))  # store trajectory information in a bundle 
# all done
t_end = time.time()
print("Time used for ", n_traj, " tracks = ", t_end - t_start)
bundle = np.array(bundle)

# save trajectory bundle in a compressed file
np.savez_compressed("bundle.npz", bundle = bundle)


#%% Plotting
# close all figures using close('all') 
# make 3d plot
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

for b in bundle:
    x = b[0]
    y = b[1]
    z = b[2]
    B.pl.plot(x,y,z, color = 'b')
B.pl.plot(xt, yt, zt, color = 'r', ls = '--')
# B.pl.plot(xpz,ypz,zpz, color = 'g')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


#%% compare 2d plot
fig2 = B.pl.figure()
ax2 = fig2.add_subplot(111)
for b in bundle:
    r = b[3]
    z = b[2]
    B.pl.plot(r,z, color = 'b')
B.pl.plot(rt,zt, color = 'r', ls = '--')
ax2.set_aspect('equal')
ax2.set_xlabel('R')
ax2.set_ylabel('Z')

#%% compare 2d plot
fig21 = B.pl.figure()
ax21 = fig21.add_subplot(111)
for b in bundle:
    x = b[0]
    y = b[1]
    B.pl.plot(x,y, color = 'b')
B.pl.plot(xt,yt, color = 'r', ls = '--')
ax21.set_aspect('equal')
ax21.set_xlabel('X')
ax21.set_ylabel('Y')

B.pl.show()

