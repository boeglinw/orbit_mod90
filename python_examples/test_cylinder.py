#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 07:52:36 2020

@author: boeglinw
"""
import calc as C
import numpy as np
import matplotlib.pyplot as pl

import BorisCylpy as BC
import acceptance_array as AA

import copy as CO

from numba import jit

# timing purposes
import timeit

import fibonacci as FI
dtr = np.pi/180.

# initialize
fg = FI.fibonacci_grid()

#%% example
# kinetic energy

K = 3*C.mev

# velocity
v = np.sqrt(2.*K/C.mp)
p_mp = C.mp*v

"""+
thv = 0.*dtr
phv = 0.*dtr

vs = np.array([np.sin(thv)*np.cos(phv),
               np.sin(thv)*np.sin(phv),
               np.cos(thv)])*v

rs = np.array([-0.0,0.,0.])
"""
Nt = 50   # number of steps in straight section

R = 0.0019
d = 0.0357
step = d/Nt


# conical detector
h_c = AA.con_detector(R, R, d)
h_c.calc_acceptance()

# largest angle
th_max = np.arctan(2.*R/d)*1.
# th_max = np.pi/2.5

# set  B-field
#Bmag = 0.3
Bmag = 0.

theta_b = 90*dtr
phi_b = 0.*dtr

# gyro radoius
#R_g = p_mp/(C.ec*Bmag)

# magnetic field direction in detector system
Bv = np.array([np.cos(phi_b)*np.sin(theta_b), np.sin(phi_b)*np.sin(theta_b), np.cos(theta_b)])*Bmag

BC.boris_cylinder.set_b0(Bv)
# initialize calc
# BC.boris_cylinder.init(1., 1836.152, R, d, step, vs )


# sunflower r theta distribution for initial position
np_pos= 20 # set the number of points
r_init, th = fg.get_circle(R, N = np_pos)
dA_pos = np.pi/np_pos * R**2

# sunflower for directions
np_dir = 200
# total solid angle
A_dir = 2.*np.pi*(1. - np.cos(th_max))
# solid angle per point
dA_dir = A_dir/np_dir
phi, theta =  fg.get_sphere(0., th_max, N = np_dir)



#%% setup result arrays and

theta_ac = []
results = []
hits = []
hits_c = []
ncalc = []


#  stop when wall is hit
BC.boris_cylinder.stop_at_hits = False
# calculate trajectory for reach position and each direction
start_time = timeit.default_timer()
for i,rr in enumerate(r_init):
    # initial position
    x = rr*np.cos(th[i])
    y = rr*np.sin(th[i])
    rs = np.array([x, y, 0.])
    for j, phv in enumerate(phi):
        thv = theta[j]
        # inital velovity
        vs = np.array([np.sin(thv)*np.cos(phv),
                       np.sin(thv)*np.sin(phv),
                       np.cos(thv)])*v
        # initialize calc
        # number of boris steps
        Nstep_l = Nt + 5
        # make sure step size and Nstep allow for a trajectory to travel the full length (d) 
        # of the detector
        # 
        step_l = step/np.cos(thv)
        BC.boris_cylinder.init(1., C.mp/C.me, R, R, d, step_l, Nstep_l, vs )
        # calc trajectory
        theta_ac.append(thv)
        nc = BC.boris_cylinder.track_cylinder(rs)
        # save results
        # results.append(CO.copy(BC.boris_cylinder.track[:,0:3]))  # save only the position
        tr = CO.copy(BC.boris_cylinder.track[:,:])
        results.append(tr)  # save position and  velocities
        hits.append(CO.copy(BC.boris_cylinder.hit))
        ncalc.append(nc)


print(f'time = {timeit.default_timer() - start_time}')


#%% make arrays

hits = np.array(hits)>0
results = np.array(results)
ncalc = np.array(ncalc)
theta_ac = np.array(theta_ac)
# hits_c = np.array(hits_c)

ok = ~hits

#%% analze results

# fraction of tracks that hit wall
ratio = np.count_nonzero(ok)/hits.shape[0]

print(f'passing fraction = {ratio:.3f}')

# calculate acceptance for the passed paricles
results_p = results[ok]
ncalc_p = ncalc[ok]

# initial and final velocities for passed particles
v_init = results_p[:,0,3:]
v_final = results_p[:,-1,3:]

v_init_mean = np.apply_along_axis(np.mean, 0, v_init)
v_final_mean = np.apply_along_axis(np.mean, 0, v_final)


# inital and final positions for passed particles
r_initial = results_p[:,0,:3]
r_final = results_p[:,-1,:3]

r_init_mean = np.apply_along_axis(np.mean, 0, r_initial)
r_final_mean = np.apply_along_axis(np.mean, 0, r_final)

# calc effective acceptance

cth_ac = np.cos(theta_ac[ok])
# cos(theta) of final direction
cth = (v_final/v)[:,2]
cthi = (v_init/v)[:,2]
  
# calculate acceptances
# using final values 
Acc = np.sum(cth)*dA_pos*dA_dir
# using initial values
Acci = np.sum(cthi)*dA_pos*dA_dir

# compare to analytical calc.
Ratio = Acc/h_c.accept

print(f'Ratio = {Acc/h_c.accept}')


"""
#%% plot inital and final positions
pl.figure()
for i,rr in enumerate(r_initial):
    r_i = rr
    r_f = r_final[i]
    pl.plot(r_i[:1], r_i[1:2], '.', color = 'g')
    pl.plot(r_f[:1], r_f[1:2],'.',  color = 'b')
pl.gca().set_aspect('equal')



#%% plot  all the results
pl.figure()
for i,res in enumerate(results):
    if hits[i] :
        pl.plot(res[:,0][:1], res[:,1][:1], '.', color = 'g')
        pl.plot(res[:,0][ncalc[i]-1:ncalc[i]], res[:,1][ncalc[i]-1:ncalc[i]],'.',  color = 'r')
    else:
        pl.plot(res[:,0][:1], res[:,1][:1], color = 'g')
        pl.plot(res[:,0][ncalc[i]-1:ncalc[i]], res[:,1][ncalc[i]-1:ncalc[i]],'.',  color = 'b')

pl.gca().set_aspect('equal')

#%% same as above but as a function of z
pl.figure()
for i,res in enumerate(results):
    if hits[i] :
        pl.plot(res[:,0][:ncalc[i]], res[:,2][:ncalc[i]], color = 'r', alpha = 0.2)
    else:
        pl.plot(res[:,0], res[:,2], color = 'b')

pl.gca().set_aspect('equal')
"""