#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
import Trpy as Tr
from mpl_toolkits.mplot3d import Axes3D

# file names should be set using these functions
#Tr.tracker.set_gfile_name('029881.00252.dat')
#track_name = 'track_11111.data'

Tr.tracker.set_gfile_name('029880.00234.dat')
track_name = 'track_21111.data'
Tr.tracker.set_efit_directory('../example_data/')

"""
Tr.limiter_control_mod.set_limiter_file_name('MASTLIMIT00.dat')
Tr.limiter_control_mod.set_limiter_directory('../example_data/')
Tr.limiter_control_mod.limiter_init()
Tr.control_mod.check_limiter=False
"""

Tr.tracker.load_flux()
track_dir = '../example_data/'


# get prev. caculated track
td = B.get_file(track_dir + track_name)
# tack locations
rt = td['r']
zt = td['z']
xt = td['x']
yt = td['y']

bphi = td['bphi']
br = td['br']
bz = td['bz']

#%%
# loop over all positions and calculate the fields
b_results = []
for i, rr  in enumerate(rt):
    zz = zt[i]
    b_results.append( Tr.em_fields_mod.bfield(rr,zz) )
b_results = np.array(b_results)

b_r = b_results[:,0]
b_z = b_results[:,1]
b_phi = b_results[:,2]

#%% retios

dev_b_z = 1. - (-b_z/bz)
dev_b_r = 1. - (-b_r/br)
dev_b_phi = 1. - (-b_phi/bphi)


print("max B_z deviation = ", np.abs(dev_b_z).max())
print("max B_r deviation = ", np.abs(dev_b_r).max())
print("max B_phi deviation = ", np.abs(dev_b_phi).max())