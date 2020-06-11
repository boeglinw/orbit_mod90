#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
import BTpy as BT
from mpl_toolkits.mplot3d import Axes3D


track_dir = '../example_data/'
track_name = 'track_11111.data'
# get prev. caculated track
td = B.get_file(track_dir + track_name)
# tack locations
rt = td['r']
zt = td['z']
xt = td['x']
yt = td['y']

#%%
# set initial valuesfor tracker
BT.tracker.particle_charge = 1.
BT.tracker.particle_mass_amu = 1.007347  
BT.tracker.particle_energy_mev = 3.
BT.tracker.track_length = 4.
BT.tracker.step_size = 0.01
BT.tracker.time_reversed = True

# file names should be set using these functions
BT.tracker.set_gfile_name('029881.00252.dat')
BT.tracker.set_efit_directory('../example_data/')
#%%
# initialize tracker: allocate space for trajectory data in BT.tracker.trajectory
BT.tracker.init_bt()
#%%
# load flux
BT.tracker.load_flux()

# set initial position
r0 = np.array([0.32623276891949599,
               1.7257147211829555,    
               2.8862357467520010E-002])

# initial velocity unit vectopr
uv0 = np.array([0.43438672503661674,
                -0.6748680841524517,      
                0.59654106489357639])
# set initial velocity
v0 = BT.tracker.vmag*uv0 


#%% setup the limiter
BT.boris_mod.set_limiter_file_name('MASTLIMIT00.dat')
BT.boris_mod.set_limiter_directory('../example_data/')
BT.boris_mod.limiter_init()
BT.boris_mod.check_limiter=True
#%%

# store results in res
res = []
for i in range(1):
    # calculate track result is in BT.tracker.trajectory
    nc = BT.tracker.get_trajectory(r0, v0)
        
    print('------ get_trajectory returned  >>', nc, ' steps calculated ')
        
    x = BT.tracker.trajectory[:nc,0]
    y = BT.tracker.trajectory[:nc,1]
    z = BT.tracker.trajectory[:nc,2]
        
    r = np.sqrt(x**2 + y**2)
        
    res.append((x, y, z,     r))  # here you could add more and more talectories and then save them as an npz file. 
# all done

#%% amek 3d plot
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

B.pl.plot(x,y,z, color = 'b')
B.pl.plot(xt, yt, zt, color = 'r', ls = '--')
# B.pl.plot(xpz,ypz,zpz, color = 'g')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


#%% compare 2d plot
fig2 = B.pl.figure()
ax2 = fig2.add_subplot(111)
B.pl.plot(r,z, color = 'b')
B.pl.plot(rt,zt, color = 'r', ls = '--')
ax2.set_aspect('equal')


B.pl.show()


