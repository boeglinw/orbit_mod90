#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:49:39 2022


Python example on how to use the is_in_limits and checklim functions

These function are in the fluxpy module

@author: boeglinw
"""

import fluxpy as FL
import numpy as np
import matplotlib.pyplot as pl

# simple boundary
ly = np.array([0., 0., 1., 1., 0.])
lx = np.array([0., 1., 1., 0., 0.])

# array of points
xp = np.array([0.5, 1.5, .2])
yp = np.array([0.1, .2, .3])
# y-axis boundary for a grid
yp2 = yp[:2]

rp = np.array([xp, yp]).T

#%% check is_in_limits

# checlk which points in  =an array of points are withing the limits lx, ly
res_lim = FL.flux.is_in_limits(rp[:,0], rp[:,1], lx, ly)

#%% check checklim 

# check which points of a grid with made of xp, yp2 combinations are inside the limits llx, ly

# checklim loops over the y-arrays for each x-array value

res_check = FL.flux.checklim(xp, yp2, lx, ly) == 1

xx,yy = np.meshgrid(xp, yp2)

# make the reslults of the same shape as the mesgh grid
is_inside = (res_check.reshape((xp.size, yp2.size))).T

# get the inside point coordinates
xx_inside = xx[is_inside]
yy_inside = yy[is_inside]


pl.plot (xx, yy, 'bo')
pl.plot (xx_inside, yy_inside, 'ro')

#%%



