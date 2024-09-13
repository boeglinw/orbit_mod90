#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:55:58 2020

Main driving script to calculate orbits or orbit bundles and prepare input files for FIDASIM

controlled by input file

this code can be run from the command line or ipython:

    %run calc_orbits.py -c 'my_control.data' -P '/where/your/modules/are/located/'

to run the example in this directory type

    %run calc_orbits.py

to get help on the arguments type

    %run calc_orbits.py --help


WB June 2024

new option: -D use detector head file that contains individual collimator geometries for each detector

default: use common geometries 

@author: boeglinw
"""


import sys
import os
import argparse as AG

import numpy as np
import copy as C

import LT.box as B   # this library can be installed using: pip install LabTools3 or cloned from GitHub: https://github.com/boeglinw/LabTools3.git
from LT.parameterfile import pfile

import pickle as PCL
import h5py

# colored trajectories
import matplotlib.colors as COL
# colors list
color_table = COL.CSS4_COLORS
color_names = list(color_table.keys())

m2cm = 100.
MeV2KeV = 1e3


# default value
reverse_velocities = False

#%% location of the python modules for orbit_mod90
orbit_mod90_python = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/python_modules'



#%% h5py helper functions

def h5py_string(s):
    return np.array(s).astype(h5py.string_dtype())

def h5py_read_string(dset):
    d_string = np.array(['']).astype(h5py.string_dtype())
    dset.read_direct(d_string)
    return d_string


#%% Save and retrive dictionaries as nos files

def save_dict(filename, dictionary, selection = None):
    """
    Save a dictionary to an npz file based on a selection 

    Parameters
    ----------
    filename : string
        filename for the output file.
    dictionary: dict
        dictionary to save
    selection : list of bools, optional
        bool arrays for the selected lines. The default is None, meaning all data.

    Returns
    -------
    None.

    """
    d = dictionary
    if selection is None:
       with open(filename, 'wb') as o:
           PCL.dump(d, o)
           o.close()
    else:
        # apply selection to all data, make temporary dictionary
        dd = dict(zip(list(d.keys()),[d[k][selection] for k in d.keys()]))
        with open(filename, 'wb') as o:
            PCL.dump(dd, o)
            o.close()


def read_dict(filename):
    f = None
    with open(filename, 'rb') as f:
        d = PCL.load(f)
        f.close()
    return d
        


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
                    default = 'data/calc_orbit_control_fidasim.data')
parser.add_argument('-P', '--orbit_mod90_path', help="Path to orbit_mod90 python modules ( will be added to PYTHONPATH",
                    default = orbit_mod90_python)

parser.add_argument('-n','--no_plot', action='store_true')

parser.add_argument('-D','--individual_detector_collimators', action='store_true')

# setup parser
args = parser.parse_args()

orbit_mod90_python = args.orbit_mod90_path

individual_detector_geometry = args.individual_detector_collimators

#%% setup system path
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


# FIDASIM data file (pickle file)

fidasim_file = cd['fidasim_file']
fidasim_name = cd['fidasim_system_name']


# make sure the name is not longer than 20 characeters
if len(fidasim_name) > 20:
    fidasim_system_name = fidasim_name[:20]
else:
    fidasim_system_name = fidasim_name



#%%  set initial values for tracker
#
Tr.tracker.particle_charge_ec = cd.get_value('particle_charge')
# Tr.tracker.particle_mass_amu = 1.007347
Tr.tracker.particle_mass_amu = cd.get_value('particle_mass_amu')

e_part_min = cd.get_value('particle_energy_mev_min')
e_part_max = cd.get_value('particle_energy_mev_max')
n_e_part = cd.get_value('n_particle_energy_steps')

e_particle = np.linspace(e_part_min, e_part_max, n_e_part)

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
trajectory_bundle_type = 'round'  # this is the only bundle type for FIDASIM at this point
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

# reverse velocities
reverse_velocities = cd.get_value('reverse_velocities')

# save inside part only
save_inside_only = cd.get_value('save_inside_only')

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
    if individual_detector_geometry:
        R_det_loc = dh['R_det'][i]; R_coll_loc = dh['R_coll'][i]; R_cyl_loc = dh['R_cyl'][i]; D_loc = dh['D'][i]
    else:
        R_det_loc = dh.par['R_det']; R_coll_loc = dh.par['R_coll']; R_cyl_loc = dh.par['R_cyl']; D_loc = dh.par['D']
        
    det_l = Det.detector(n,
                        dh['Detector_name'][i],
                        head_position = np.array([R_p, Z_p, Phi_p]),
                        det_pos_local = np.array([dh['xd'][i], dh['yd'][i], dh['zd'][i]]),
                        det_direction = np.array([dh['theta_d'][i]*dtr, dh['phi_d'][i]*dtr]),
                        arm_rotation = arm_rotation,
                        tracker = Tr,
                        bundle_fname = f'det_{n}_'+dh['Detector_name'][i]+'.npz',
                        color = dh['color'][i],
                        R_det = R_det_loc, R_coll = R_coll_loc, R_cyl = R_cyl_loc, D = D_loc,
                        zero_at_coll = zero_at_collimator,
                        fib_scale = fib_angle_scale,
                        ignore_detector_field = ignore_detector_bfield)
    detector_head.append(det_l)

# calculate trajectories

for det_l in detector_head:
    det_l.all_central = []
    det_l.all_bundles = []
    det_l.all_Bf_bundles = []
    det_l.all_hdf_info = []
    det_l.all_acc = []
    det_l.get_fidasim_geometry()
    for E_p in e_particle:
        Tr.tracker.particle_energy_mev = E_p  # set new particle energy in MeV
        Tr.tracker.init_tracker()             # init tracker wirh new energy          
        det_l.init_trajectories(N_pos = N_pos_det, N_dir = N_dir_det)   # init trajectory calculation
        det_l.calc_trajectories(bundle_type = trajectory_bundle_type,save = False, inside_only = save_inside_only)
        det_l.calc_central(save = False)
        # calculations complete store results
        det_l.setup_hdf_data()
        # store results
        det_l.all_hdf_info.append([det_l.hdf_nsteps, det_l.hdf_n_steps_actual,det_l.hdf_n_rays])
        det_l.all_central.append(C.copy(det_l.central_track))  # store central tracks
        det_l.all_bundles.append(C.copy(det_l.bundle))         # store trajectory bundles (rays, sightlines)
        det_l.all_Bf_bundles.append(C.copy(det_l.Bf_bundle))   # store corresponding B-fields 
        det_l.all_acc.append(C.copy(det_l.acc))
if args.no_plot:
    sys.exit()
    
    
    
#%% prepare data for HDF input structure
# create a system name

# array of actual step numbers

na = np.array([[hi[1] for hi in det_l.all_hdf_info] for det_l in detector_head])
# arrange the axis so that the sequence of indices is as is read in FIDASIM

nactual_a = np.moveaxis(na, 0, 2)

#%% array of acceptances
# array of actual step numbers
daomega = np.array([[acc for acc in det_l.all_acc] for det_l in detector_head])*m2cm**2

# arrange the axis so that the sequence of indices is as is read in FIDASIM
daomega_a = np.moveaxis(daomega, 0, -1)


#%% sightline array

# array of indices for selected trajectory variables
p_array = np.array([detector_head[0].bundle_vars[k] for k in ['vr', 'vphi', 'vz', 'r','phi','z' ]])
pb_array = np.array([detector_head[0].Bf_bundle_vars[k] for k in ['b_pol_r', 'b_phi', 'b_pol_z', 'b_total','psi_rel']])

# indices into traj_data
v_array = np.array([0,1,2])
pos_var_array = np.array([0,1,2,3,5])  # vr, vphi, vz, r, z
dsightline = []


s_data = []
sb_data = []
# loop over all detector heads
for det_l in detector_head[:]:
    e_data = []
    eb_data = []  # B-field data array
    # loop over all detector trajectory bundles
    for i,e_bundle in enumerate(det_l.all_bundles[:]):
        t_data = []
        tb_data = []
        e_Bf_bundle = det_l.all_Bf_bundles[i]
        # loop over all different energies 
        for j,t in enumerate(e_bundle[:]):
            bf_loc = e_Bf_bundle[j]
            ns = det_l.all_hdf_info[i][1][j]  # get the corresponding number of steps
            # copy trajectory to traj_data
            traj_data = np.zeros(shape = (Tr.tracker.nsteps, 6))
            trbf_data = np.zeros(shape = (Tr.tracker.nsteps, 5))
            traj_data[0:ns,:] = t[:,p_array]
            trbf_data[0:ns,:] = bf_loc[:,pb_array]
            # reverse velocity direction if selected
            if reverse_velocities:
                traj_data[0:ns,v_array] *= -1
            # convert position dependent variables to cm
            traj_data[0:ns,pos_var_array] *= m2cm
            t_data.append(traj_data)
            tb_data.append(trbf_data)
        e_data.append(t_data)
        eb_data.append(tb_data)
    s_data.append(e_data)
    sb_data.append(eb_data)
    
ds_data = np.array(s_data)
dsb_data = np.array(sb_data)

# move axes to be compatible with FIDASIM input  
sightline_a = np.moveaxis(ds_data, [0,1,-1,-2], [-1,0,1,2])

#%% create dictionary for FIDASIM input

FD = {}

nchan = len(detector_head)

FD['system'] = fidasim_system_name  # max 20 characters
FD['id'] = [det_l.name for det_l in detector_head ]
FD['d_shape'] = nchan*['round']
FD['a_shape'] = nchan*['round']
FD['a_cent'] = np.array([dd.a_c for dd in detector_head]).T * m2cm
FD['a_redge'] = np.array([dd.a_r for dd in detector_head]).T * m2cm
FD['a_tedge'] = np.array([dd.a_t for dd in detector_head]).T * m2cm
FD['d_cent'] = np.array([dd.d_c for dd in detector_head]).T * m2cm
FD['d_redge'] = np.array([dd.d_r for dd in detector_head]).T * m2cm
FD['d_tedge'] = np.array([dd.d_t for dd in detector_head]).T * m2cm 
FD['nchan'] = nchan
FD['nenergy'] = e_particle.shape[0]
FD['nrays'] = det_l.bundle.shape[0] 
FD['nsteps'] = Tr.tracker.nsteps + 0
FD['earray'] = e_particle*MeV2KeV 
FD['nactual'] = nactual_a
FD['daomega'] = daomega_a
FD['sightline'] = sightline_a
FD['sightline_Bf'] = dsb_data

# save dictionary
save_dict(fidasim_file, FD)

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
def plot_plasma3D(central = False, plot_all = False):
    # close all figures using close('all')
    # make 3d plot
    fig3d = B.pl.figure()
    ax = fig3d.add_subplot(111, projection='3d')
    ax.plot_surface(XX, YY, ZZ, color = 'r', alpha = 0.2)
    # plot all  in det_l
    # plot all  in detector_head:
    if plot_all:
        for dd in detector_head:
            if central:
                for ddd in dd.all_central:
                    x = ddd[:,0]
                    y = ddd[:,1]
                    z = ddd[:,2]
                    B.pl.plot(x,y,z, color = color_table[dd.color])
            else:
                for ddd in dd.all_bundles:
                    for i,b in enumerate(ddd):
                        x = b.T[0]
                        y = b.T[1]
                        z = b.T[2]
                        B.pl.plot(x,y,z, color = color_table[dd.color])
    else:
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

    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


#%% RZ plot
def plot_plasma_rz(central = False, plot_all = False):
    fig_rz = B.pl.figure()
    limiter.draw_side_all()
    
    # draw the flux
    B.pl.contour(rrg,zzg, psi, levels = 40)
    
    # draw the plasma
    B.pl.plot(rbdry, zbdry, color = 'r')

    # plot all  in detector_head:
    if plot_all:
        for dd in detector_head:
            if central:
                for ddd in dd.all_central:
                    r = ddd[:,3]
                    z = ddd[:,2]
                    B.pl.plot(r,z, color = color_table[dd.color])
            else:
                for ddd in dd.all_bundles:
                    for i,b in enumerate(ddd):
                        r = b.T[3]
                        z = b.T[2]
                        B.pl.plot(r,z, color = color_table[dd.color])
    else:
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
def plot_plasma_mid(central = False, plot_all = False):
    fig_mid = B.pl.figure()
    limiter.draw_top_all()
    
    # draw the plasma
    plot_ring(r_plasma_min, r_plasma_max, B.pl.gca(), color = 'r', alpha = 0.2, edgecolor = None)
    
    # plot all  in detector_head:
    if plot_all:
        for dd in detector_head:
            if central:
                for ddd in dd.all_central:
                    x = ddd[:,0]
                    y = ddd[:,1]
                    B.pl.plot(x,y, color = color_table[dd.color])
            else:
                for ddd in dd.all_bundles:
                    for i,b in enumerate(ddd):
                        x = b.T[0]
                        y = b.T[1]
                        B.pl.plot(x,y, color = color_table[dd.color])
    else:
        for dd in detector_head:
            if central:
                x = dd.central_track.T[0]
                y = dd.central_track.T[1]
                B.pl.plot(x,y, color = color_table[dd.color])
            else:            
                for i,b in enumerate(dd.bundle):
                    x = b.T[0]
                    y = b.T[1]
                    B.pl.plot(x,y, color = color_table[dd.color])    

    B.pl.ylim((-2.,2.))
    B.pl.xlim((-2.,2.))
    
