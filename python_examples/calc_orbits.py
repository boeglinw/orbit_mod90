#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

Main driving script to calculate orbits or orbit bundles

controlled by input file

this code can be run from the command line or ipython:

    %run calc_orbits.py -c 'my_control.data' -P '/where/your/modules/are/located/'

to run the example in this directory type

    %run calc_orbits.py

to get help on the arguments type

    %run calc_orbits.py --help

@author: boeglinw
"""


import sys
import os
import argparse as AG

import numpy as np

import LT.box as B   # this library can be installed using: pip install LabTools3 or cloned from GitHub: https://github.com/boeglinw/LabTools3.git
from LT.parameterfile import pfile

# colored trajectories
import matplotlib.colors as COL
# colors list
color_table = COL.CSS4_COLORS
color_names = list(color_table.keys())




#%% setup the PYTHON path

def add_sys_path(new_path):
    """ AddSysPath(new_path): adds a directory to Python's sys.path

    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os

    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1

    # Standardize the path. Windows is case-insensitive, so lowercase
    # for definiteness.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower(  )

    # Check against all currently available paths
    for x in sys.path:
        x = os.path.abspath(x)
        if sys.platform == 'win32':
            x = x.lower(  )
        if new_path in (x, x + os.sep):
            return 0
    sys.path.append(new_path)
    return 1


#%% plot a ring, used for plasma plotting
def plot_ring(rs, rl, ax, n = 50 , **kwargs ):
    radii = [rs, rl]
    theta = np.linspace(0, 2*np.pi, n, endpoint=True)
    xs = np.outer(radii, np.cos(theta))
    ys = np.outer(radii, np.sin(theta))

    # in order to have a closed area, the circles
    # should be traversed in opposite directions
    xs[1,:] = xs[1,::-1]
    ys[1,:] = ys[1,::-1]

    ax.fill(np.ravel(xs), np.ravel(ys), **kwargs)# , edgecolor='#348ABD')



#%% input section get comman line arguments
parser = AG.ArgumentParser(formatter_class=AG.ArgumentDefaultsHelpFormatter)

HERE = os.getcwd()

parser.add_argument('-c', '--control_file', help="orbit calculation control file",
                    default = 'calc_orbit_control.data')
parser.add_argument('-P', '--orbit_mod90_path', help="Path to orbit_mod90 python modules ( will be added to PYTHONPATH",
                    default = HERE+'/../python_modules')

parser.add_argument('-n','--no_plot', action='store_true')

# setup parser
args = parser.parse_args()

orbit_mod90_python = args.orbit_mod90_path

#%% location of the python modules for orbit_mod90
orbit_mod90_python = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/python_modules'

if (add_sys_path(orbit_mod90_python)) < 0 :
    print(f'Cannot add {orbit_mod90_python} to PYTHONPATH as it does not exist !')
    sys.exit(-1)

# setup orbit specific modules
# Tracking
import Trpy as Tr

# detectors
import detectors as Det

# limiter drawing
import get_limiter as gl

dtr = np.pi/180.

trackers = {'Boris':Tr.tracker.boris_t, 'Bulirsch_stoer':Tr.tracker.bulirsch_stoer_t}


#%%  control file input

# load control file, this is a pure parameter file
cd = pfile(args.control_file)

# this is a data file describing the detector head
dh = B.get_file(cd['detector_head'])


#%%  set initial values for tracker
#
Tr.tracker.particle_charge_ec = cd.get_value('particle_charge')
# Tr.tracker.particle_mass_amu = 1.007347
Tr.tracker.particle_mass_amu = cd.get_value('particle_mass_amu')
Tr.tracker.particle_energy_mev = cd.get_value('particle_energy_mev')

Tr.tracker.track_length = cd.get_value('track_length')
Tr.tracker.step_size = cd.get_value('step_size')
Tr.control_mod.time_reversed = cd.get_value('time_reversed')
Tr.control_mod.reverse_poloidal_flux = cd.get_value('reverse_poloidal_flux')

# Select which tracker to use
Tr.tracker.selected_tracker = trackers[cd.get_value('selected_tracker')]

# file names should be set using these functions
Tr.tracker.set_gfile_name(cd.get_value('gfile_name'))
Tr.tracker.set_efit_directory(cd.get_value('efit_directory'))


#%% Control trajectory bundle calculation
trajectory_bundle_type = cd.get_value('bundle_type')
# number of fibonacci positions in detector area
N_pos_det = cd.get_value('number_of_detector_points')
# number of fibonacci points for direction
N_dir_det = cd.get_value('number_of_directions')
# set collimaor at z = 0 for track
zero_at_collimator = cd.get_value('zero_at_collimator')
# scale factor for fibonacci angle grid
fib_angle_scale = cd.get_value('fib_angle_scale')

# ingore b-field
ignore_detector_bfield = cd.get_value('ignore_detector_bfield')

#%% initialize tracker: allocate space for trajectory data in BT.tracker.trajectory
#
Tr.tracker.init_tracker()
#%% load flux i.e the magnetic field and its interpolations
#
Tr.tracker.load_flux()


#%% setup the limiter
Tr.limiter_control_mod.set_limiter_file_name(cd.get_value('limiter_file_name'))
Tr.limiter_control_mod.set_limiter_directory(cd.get_value('limiter_directory'))
Tr.limiter_control_mod.limiter_init()
Tr.control_mod.check_limiter = cd.get_value('check_limiter')
Tr.control_mod.print_hit = cd.get_value('print_hit')
# setup for limiter drawing
limiter = gl.limiter('limiter_drawing.data')
# Tr.control_mod.print_polygon = False

#%% setup a detectors
R_p = cd.get_value('detector_head_R')
Z_p = cd.get_value('detector_head_Z')
Phi_p = cd.get_value('detector_head_Phi')*dtr


arm_rotation = cd.get_value('arm_rotation')*dtr

detector_head = []

for i,n  in enumerate(dh['Detector_number']):
    det_l = Det.detector(n,
                        dh['Detector_name'][i],
                        head_position = np.array([R_p, Z_p, Phi_p]),
                        det_pos_local = np.array([dh['xd'][i], dh['yd'][i], dh['zd'][i]]),
                        det_direction = np.array([dh['theta_d'][i]*dtr, dh['phi_d'][i]*dtr]),
                        arm_rotation = arm_rotation,
                        tracker = Tr,
                        bundle_fname = f'det_{n}_'+dh['Detector_name'][i]+'.npz',
                        color = dh['color'][i],
                        R_det = dh.par['R_det'], R_coll = dh.par['R_coll'], R_cyl = dh.par['R_cyl'], D = dh.par['D'],
                        zero_at_coll = zero_at_collimator,
                        fib_scale = fib_angle_scale,
                        ignore_detector_field = ignore_detector_bfield)
    detector_head.append(det_l)
# initialize
for det_l in detector_head:
    det_l.init_trajectories(N_pos = N_pos_det, N_dir = N_dir_det)

# calculate trajectories
for det_l in detector_head:
    det_l.calc_trajectories(bundle_type = trajectory_bundle_type)
    det_l.calc_central()

if args.no_plot:
    sys.exit()
#%% prepare to plot the plasma
# get the boundary data
nbdry = Tr.flux_par_mod.nbdry
zbdry = Tr.flux_par_mod.zbdry[:nbdry]
rbdry = Tr.flux_par_mod.rbdry[:nbdry]

phibdry = np.linspace(0., 2.*np.pi, 36)
RR,PP = np.meshgrid(rbdry, phibdry)
ZZ,PP = np.meshgrid(zbdry, phibdry)
XX = RR*np.cos(PP)
YY = RR*np.sin(PP)

r_plasma_min = rbdry.min()
r_plasma_max = rbdry.max()

#%% prepare plotting the flux
mh = Tr.flux_par_mod.mh # number of horzontal flux grid point
mw = Tr.flux_par_mod.mw # number of vertical flux grid points
n_tot_flux = mw*mh      # number of total flux grid points

# get the flux as 2d array
psi = Tr.flux_par_mod.psi[:n_tot_flux].reshape((mw,mh)).T

# setup the grid point arrays
rg = Tr.flux_par_mod.rgrid[:mw]
zg = Tr.flux_par_mod.zgrid[:mw]

rrg,zzg = np.meshgrid(rg, zg)


#%% Plotting
def plot_plasma3D(central = False):
    # close all figures using close('all')
    # make 3d plot
    fig3d = B.pl.figure()
    ax = fig3d.add_subplot(111, projection='3d')
    ax.plot_surface(XX, YY, ZZ, color = 'r', alpha = 0.2)
    # plot all  in det_l
    for dd in detector_head:
        if central:
            x = dd.central_track.T[0]
            y = dd.central_track.T[1]
            z = dd.central_track.T[2]
            B.pl.plot(x,y,z, color = color_table[dd.color])
        else:
            for b in dd.bundle:
                x = b.T[0]
                y = b.T[1]
                z = b.T[2]
                B.pl.plot(x,y,z, color = color_table[dd.color])
    # B.pl.plot(xpz,ypz,zpz, color = 'g')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


#%% RZ plot
def plot_plasma_rz(central = False):
    fig_rz = B.pl.figure()
    limiter.draw_side_all()
    
    # draw the flux
    B.pl.contour(rrg,zzg, psi, levels = 40)
    
    # draw the plasma
    B.pl.plot(rbdry, zbdry, color = 'r')

    # plot all  in detector_head:
    for dd in detector_head:
        if central:
            r = dd.central_track.T[3]
            z = dd.central_track.T[2]
            B.pl.plot(r,z, color = color_table[dd.color])
        else:            
            for i,b in enumerate(dd.bundle):
                r = b.T[3]
                z = b.T[2]
                B.pl.plot(r,z, color = color_table[dd.color])
    B.pl.ylim((-2.,2.))
    B.pl.xlim((0.,2.))



#%% midplane plot
def plot_plasma_mid(central = False):
    fig_mid = B.pl.figure()
    limiter.draw_top_all()
    
    # draw the plasma
    plot_ring(r_plasma_min, r_plasma_max, B.pl.gca(), color = 'r', alpha = 0.2, edgecolor = None)
    
    # plot all  in det_l
    for dd in detector_head:
        if central:
            x = dd.central_track.T[0]
            y = dd.central_track.T[1]
            B.pl.plot(x,y, color = color_table[dd.color])
        else:               
            for b in dd.bundle:
                x = b.T[0]
                y = b.T[1]
                B.pl.plot(x,y, color = color_table[dd.color])
    B.pl.ylim((-2.,2.))
    B.pl.xlim((-2.,2.))
    
