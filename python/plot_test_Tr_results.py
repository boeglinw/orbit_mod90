#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:44:24 2020

@author: boeglinw
"""

import numpy as np
import LT.box as B
from mpl_toolkits.mplot3d import Axes3D

def rotate(r, phi):
    rx = r[0]*np.cos(phi) - r[1]*np.sin(phi)
    ry = r[0]*np.sin(phi) + r[1]*np.cos(phi)
    return (rx,ry)

    

track_dir = '../example_data/'
#track_name = 'track_11111.data'
track_name = 'test_run_track_11111.data'


# get prev. caculated track
td = B.get_file(track_dir + track_name)
# tack locations
rt = td['r']
zt = td['z']
xt = td['x']
yt = td['y']

# fields
brt = td['br']
bzt = td['bz']
bphit = td['bphi']



# boris 
d0 = np.loadtxt('track_g_Tr_BO.dat')
xb = d0[:,0]
yb = d0[:,1]
zb = d0[:,2]
rb = np.sqrt(xb**2 + yb**2)


# bulirsch calc.
d = np.loadtxt('track_g_Tr_BS_cart.dat')
#d = np.loadtxt('track_g_3.0.dat')
x = d[:,0]
y = d[:,1]
z = d[:,2]
r = np.sqrt(x**2 + y**2)

fig2 = B.pl.figure()
ax2 = fig2.add_subplot(111)
B.pl.plot(r,z, color = 'b')
B.pl.plot(rb,zb, color = 'g')
B.pl.plot(rt,zt, color = 'r', ls = '--')
ax2.set_aspect('equal')

fig3 = B.pl.figure()
ax3 = fig3.add_subplot(111)
B.pl.plot(x, y, color = 'b')
B.pl.plot(xb,yb, color = 'g')
B.pl.plot(xt,yt, color = 'r', ls = '--')
ax3.set_aspect('equal')

"""
dpz = np.loadtxt('track_g_3.0_pz.dat')
xpz = dpz[:,0]
ypz = dpz[:,1]
zpz = dpz[:,2]
rpz = np.sqrt(xpz**2 + ypz**2)
"""

#%%
fig3d = B.pl.figure()
ax = fig3d.add_subplot(111, projection='3d')

B.pl.plot(x,y,z, color = 'b')
B.pl.plot(xb,yb,zb, color = 'g')
B.pl.plot(xt, yt, zt, color = 'r', ls = '--')
# B.pl.plot(xpz,ypz,zpz, color = 'g')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

    