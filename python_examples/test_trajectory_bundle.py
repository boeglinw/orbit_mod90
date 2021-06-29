#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 07:52:36 2020

@author: boeglinw
"""
import calc as C
import numpy as np
import matplotlib.pyplot as pl


# detectors
import detectors as Det
import acceptance_array as AA

import copy as CO


dtr = np.pi/180.

#%% example
# kinetic energy

K = 3*C.mev

# velocity
v = np.sqrt(2.*K/C.mp)

mp_me = C.mp/C.me

"""
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


# sunflower r theta distribution for initial position
np_pos= 50 # set the number of points

# sunflower for directions

th_scale = 3
np_dir = int(100*th_scale**2)


bd = MTb.bundle(  N_pos = np_pos,  # number of fibonacci points in detector area
                         N_dir = np_dir, # number of direction fibonacci points
                         N_s = Nt, dN_s = 5,        # number of steps inside detector
                         v = v,                    # particle velocity
                         m_me = mp_me,
                         q_ec = 1.,
                         R_det = Rd, R_coll = Rc, R_cyl = Rcyl, D = d, # detector geometry
                         Bf = Bv,                           #  mag. field at detector location
                         scale = th_scale,
                         z_zero = False)



#%% calculations
bd.initialize()
bd.calc_bundle()
rf,vf,A = bd.get_bundle()
r_f = rf.T
v_f = vf.T

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



