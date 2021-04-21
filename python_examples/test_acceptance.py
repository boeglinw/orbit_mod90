#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:08:24 2020


Test acceptance calculation s

@author: boeglinw
"""

import numpy as np

# acceptance module
import acceptance_array as AA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import calc as C

import fibonacci as FI

from numba import jit

# timing purposes
import timeit

import LT.box as B

twopi = 2.*np.pi
dtr = np.pi/180.

def get_angle(x,y):
    th = np.arctan2(y,x)
    sel = th < 0.
    th[sel] = th[sel] + twopi
    return th



def plot_vector(r):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    r0 = (0., 0., 0.)
    ax.quiver(*r0, *r)
    ax.set_xlim([-1, 1.])
    ax.set_ylim([-1, 1.])
    ax.set_zlim([-1, 1.])
def plot_points(x,y,z, xlim = [-1.,1.], ylim = [-1.,1.], zlim = [-1., 1.]):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, '.')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

def plot_lines(x,y,z, xlim = [-1.,1.], ylim = [-1.,1.], zlim = [-1., 1.]):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

# setup a collimator

# hole detector
Rh = .002
D = .05

Rc = 0.002
Rd = 0.003

Ah = np.pi*Rh**2

# rectangular collimator
C_s = np.sqrt(Ah)/2.  # coll. half width
D_s = np.sqrt(Ah)/2.  # det. half width


# square detector
r_D = AA.rect_detector(C_s, C_s, D_s, D_s, D)
r_D.calc_acceptance()



h_D = AA.cyl_detector(Rh, D)
h_D.calc_acceptance()

# acc_geo = (np.pi*Rh**2)**2/ (4.*np.pi*D**2)
acc_geo = (np.pi*Rh**2)**2/ (D**2)

A = np.pi*Rh**2


# conical detector
h_c = AA.con_detector(Rh, Rh, D)
h_c.calc_acceptance()

h_c1 = AA.con_detector(Rc, Rd, D)
h_c1.calc_acceptance()



#%%
# calculate acceptance using MC
#
N_events = 1000
# largest angle
th_max_mc = np.arctan(2.*Rh/D)
cth_max_mc = np.cos(th_max_mc)
# select random values for cos(th) to create uniformly distributed points
cth_mc = np.random.uniform(low = cth_max_mc, high = 1., size = N_events)
th_mc = np.arccos(cth_mc)
# select phi randomlu
phi_mc = np.random.uniform(low = 0., high = 2.*np.pi, size = N_events)
# Full solid angle
delta_omega_mc = 2.*np.pi*(1. - cth_max_mc)
# calculate overlap
Ov_mc = AA.H_overlap(th_mc, Rh, D)*np.cos(th_mc)
# acceptance
Ac = np.sum(Ov_mc)*delta_omega_mc/N_events

print(f'ratio = {Ac/h_D.accept}')

print (f'geometrical acceptance ratio = ')

#%% plot
xp = np.cos(phi_mc)*np.sin(th_mc)
yp = np.sin(phi_mc)*np.sin(th_mc)
zp = np.cos(th_mc)

xy_lim = np.array([xp.min(), xp.max()])
z_lim = np.array(zp.min(), zp.max())
plot_points(xp, yp, zp, xlim = xy_lim, ylim = xy_lim, zlim = z_lim)

#%% test sub divisition
# create a conical detector
#hc = AA.con_detector(Rh, Rh, D)
hc = AA.con_detector(Rc, Rd, D)
hc.calc_acceptance()

# of cos theta steps
nth = 50

# number of phi steps
nph = 20

# create an array if equally spaced cos(theta) values
ctha = np.linspace(np.cos(hc.th_max*1.1), 1., nth)[::-1]
tha = np.arccos(ctha)
# tha = np.linspace(0., hc.th_max, nth)
dth = np.diff(tha)
# dth = np.pad(dth, (0,1), 'edge')  # repeat last element

# array of phi values
ph = np.linspace(0., 2.*np.pi, nph)
dph = np.diff(ph)
# dph = np.pad(dph, (0,1), 'edge')  # repeat last element

# create all possible theta phi combinations
P,T = np.meshgrid(ph, tha[:-1])



#%% integrate using simpsons rule
# flatten arrays
Tf = T.flatten()
Pf = P.flatten()

# plot the points
#pl.polar(Pf, Tf, '.')

#%% plot
xp = np.cos(Pf)*np.sin(Tf)
yp = np.sin(Pf)*np.sin(Tf)
zp = np.cos(Tf)

xy_lim = np.array([xp.min(), xp.max()])
z_lim = np.array(zp.min(), zp.max())
# z_lim = 1. + xy_lim

plot_points(xp, yp, zp, xlim = xy_lim, ylim = xy_lim, zlim = z_lim)


#%%  calculate the overlap function

# 
fcalc = hc.f(Tf)*Rc**2
# reshape array
fcalca = fcalc.reshape(T.shape) 
# integrate for constant phi using simpson
theta_integral = AA.integrate.simps(fcalca, T, axis = 0)

# print(f'theta ratio = {theta_integral/theta_sum}')
accept_t = AA.integrate.simps(theta_integral, ph)



print(f'ratio sub = {accept_t/hc.accept}')


#%% uniformly distributed points Fibonacci Grid

num_pts = 2000
fg = FI.fibonacci_grid(num_pts)

#indices = np.arange(0, num_pts, dtype=float) + 0.5
#arg = 1 - 2*indices/num_pts

# largest angle
#th_max = np.arctan((Rd + Rc)/D)
th_max = np.arctan((Rh + Rh)/D)

# total solid angle
A = 2.*np.pi*(1. - np.cos(th_max))
# solid angle per point
dA = A/num_pts

phi, theta = fg.get_sphere(0., th_max)


xp = np.cos(phi)*np.sin(theta)
yp = np.sin(phi)*np.sin(theta)
zp = np.cos(theta)


# over = AA.C_overlap(theta, Rc, Rd, D)
over = AA.C_overlap(theta, Rh, Rh, D)

acc_a = over*np.cos(theta)*dA

accept_f = acc_a.sum()

print(f'ratio = {accept_f/h_c.accept}')

#%% sunflower
# simple version: see https://demonstrations.wolfram.com/SunflowerSeedArrangements/
#

Nr = 50
r, th = fg.get_circle(Rh, N = Nr)
dA_r = np.pi*Rh**2/Nr

# modified version from https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle

#%% straight trajectory
@jit(nopython=True)
def calc_straight(r, th, ph , Ns, d):
    dz = d/Ns
    uv = np.array([np.cos(ph)*np.sin(th), np.sin(ph)*np.sin(th), np.cos(th)])
    args = np.arange(Ns + 1)
    z = args*dz
    l = z/np.cos(th)
    ls = l.reshape((l.shape[0],1))
    uvs = uv.reshape((1, uv.shape[0]))
    t = np.dot(ls, uvs) + r
    return t


#%%  straight trajectories
Nt = 50   # number of steps in straight section
tracks = []
hits = []
theta_t = []
# loop over all positions and directions
for i, r_r in enumerate(r):
    th_r = th[i] 
    r0 = np.array([r_r*np.cos(th_r), r_r*np.sin(th_r), 0.])
    for j, theta_s in enumerate(theta): 
        phi_s = phi[j]
        tr = calc_straight(r0, theta_s, phi_s, Nt, D)
        rt = np.sqrt(tr[:,0]**2 + tr[:,1]**2)
        hit = (rt>Rh).max()
        tracks.append(tr)
        hits.append(hit)
        theta_t.append(theta_s)
        
tracks = np.array(tracks)
hits = np.array(hits)
theta_t = np.array(theta_t)
# estimate solid angle

cth = np.cos(theta_t[~hits])

Acc = np.sum(cth)*dA_r*dA


Ratio = Acc/h_c.accept

print(f'Ratio = {Acc/h_c.accept}')

#%% calculate trajectory in constant B- field in y-direction
@jit(nopython=True)
def Ef(x):
    return np.array([0., 0., 0.])

@jit(nopython=True)
def Bf(x):
    return np.array([0., 0., 0.])

@jit(nopython=True)
def calc_trajectory(x0, theta, phi, q = C.ec, m = C.mp, E=3.*C.mev, l=0.1, dl = 1e-2):
    # calculate trajectroy between detector and collimator
    # initial velocity unit vector
    u = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
    # particle properties
    qoverm = q/m     
    vmag = np.sqrt(2.*E/m)
    V0 = vmag*u
    # trajectory length and step size
    dt = dl/vmag
    Nsteps = int(l/dl)
    xreturn = np.zeros(shape = (Nsteps+1, 3))
    Vreturn = np.zeros(shape = (Nsteps+1, 3))
    # EM fields
    for i in range(Nsteps+1):
        Vminus = V0 + qoverm*dt/2.*Ef(x0)
        tvect = qoverm*dt/2.*Bf(x0)
        Vprime = Vminus + np.cross(Vminus,tvect)
        svect = 2.*tvect/(1. + (tvect**2).sum())
        Vplus = Vminus + np.cross(Vprime, svect)
        V1 = Vplus + qoverm*dt/2.*Ef(x0)
        x1 = x0 + V1*dt
        xreturn[i,:] = x0
        Vreturn[i,:] = V0
        x0 = x1
        V0 = V1
    return xreturn, Vreturn

#%% calculate trajectory in constant B- field in y-direction in a normalized fashiom

@jit(nopython=True)
def calc_trajectory_n(x0, theta, phi, q = C.ec, m = C.mp, E=3.*C.mev, l=0.1, dl = 1e-2):
    # calculate trajectroy between detector and collimator
    # initial velocity unit vector
    u = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
    # particle properties
    qoverm = q/m     
    vmag = np.sqrt(2.*E/m)
    V0 = vmag*u
    # trajectory length and step size
    dt = dl/vmag
    qp = qoverm*dt/2.
    Nsteps = int(l/dl)
    xreturn = np.zeros(shape = (Nsteps+1, 3))
    Vreturn = np.zeros(shape = (Nsteps+1, 3))
    # EM fields
    for i in range(Nsteps+1):
        Vminus = V0 + qp*Ef(x0)
        tvect = qp*Bf(x0)
        Vprime = Vminus + np.cross(Vminus,tvect)
        svect = 2.*tvect/(1. + (tvect**2).sum())
        Vplus = Vminus + np.cross(Vprime, svect)
        V1 = Vplus + qp*Ef(x0)
        x1 = x0 + V1*dt
        xreturn[i,:] = x0
        Vreturn[i,:] = V0
        x0 = x1
        V0 = V1
    return xreturn, Vreturn


#%% calc curved trajectory
class kinematics:
    def __init__(self,  q = C.ec, m = C.mp, E=3.*C.mev, l=0.1, dl = 1e-2, Bf = np.array([0.,.0, 0.])):
                 
                 self.q = q
                 self.m = m
                 self.E = E
                 self.l = l
                 self.dl = dl
                 self.B = Bf
                 
# thats it
setup0 =  kinematics( E = 3.*C.mev, l = 1.0*D, dl = D/50., Bf = np.array([0.,.0, 0.]))

#%%  calc tracks
tracks_c = []
vtracks_c = []
hits_c = []
theta_t_c = []
start_time = timeit.default_timer()
# loop over all positions and directions
for i, r_r in enumerate(r):
    th_r = th[i] 
    r0 = np.array([r_r*np.cos(th_r), r_r*np.sin(th_r), 0.])
    for j, theta_s in enumerate(theta): 
        phi_s = phi[j]
        r_tr, v_tr = calc_trajectory_n(r0, theta_s, phi_s, E = 3.*C.mev, l = 1.0*D, dl = D/50. )
        rt = np.sqrt(r_tr[:,0]**2 + r_tr[:,1]**2)
        hit = (rt>Rh).max()
        tracks_c.append(r_tr)
        vtracks_c.append(v_tr)
        hits_c.append(hit)
        theta_t_c.append(theta_s)
        
tracks_c = np.array(tracks_c)
hits_c = np.array(hits_c)
theta_t_c = np.array(theta_t_c)
print(f'time = {timeit.default_timer() - start_time}')
# estimate solid angle

cth_c = np.cos(theta_t_c[~hits_c])

Acc_c = np.sum(cth_c)*dA_r*dA


Ratio = Acc_c/h_c.accept

print(f'Ratio = {Acc_c/h_c.accept}')


    