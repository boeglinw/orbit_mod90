#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:28:36 2020

@author: boeglinw
"""
import sympy as S

# setup rotation matrices
# rotation around x-axis by x
Rx = lambda x: S.Matrix([[1,0,0],[0, S.cos(-x), -S.sin(-x)], [0, S.sin(-x), S.cos(-x)]])


# rotation around y-axis by x
Ry = lambda x: S.Matrix([[S.cos(x),0,-S.sin(x)],[0, 1, 0], [S.sin(x), 0 , S.cos(x)]])

# rotation around z-acis by x
Rz = lambda x: S.Matrix([[S.cos(x),-S.sin(x),0], [S.sin(x), S.cos(x), 0], [0, 0, 1]])


# Total transformation from local detector system to local probe stsrem
def T(t,p):
    return Ry(t)*Rx(p)



# velocity vector
vx, vy, vz = S.symbols('v_x v_y v_z')
v = S.Matrix([[vx, vy, vz]]).T
#%%
# calculate rotation matrices for special cases
# examples

vt = S.Matrix([0,0,1])   # unitvector along z-axis
# example1
theta = S.pi/2
phi = 0
T(theta, phi)*vt

#%%
theta = S.pi/2
phi = S.pi/4
T(theta, phi)*vt

#%%
theta = S.pi/4
phi = S.pi/4
T(theta, phi)*vt


# transformation to overall coord. system
theta = S.pi/4
phi = S.pi/4
vl = T(theta, phi)*vt # vector in local coordinate sysem

v_g = Rz(S.pi/2)*vl
