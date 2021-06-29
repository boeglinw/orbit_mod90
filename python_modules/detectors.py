#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 07:11:37 2021

Detector module describing each detector, location and orientation. Is is also used to calculate trajectories

Each detector is part of a detector head. The so-called local coordinate system is identical with the detector head coordinate system. 
Detector positions are most frequently given with respect to the local (head) coordinate system. Each detector has is own local coordinate system that is
used to generate tracks (see coordinate_systems module)

@author: boeglinw
"""
import copy as C
import time
import numpy as np


# coordinate transformations
import coordinate_systems as Cs
# detector-collimator trajectory setup
import make_trajectory_bundles as MTb


dtr = np.pi/180.


#%% detector class

class detector:
    def __init__(self, number,  name,
                 R_det = 0.001,
                 R_coll = 0.001,
                 R_cyl = 0.001,
                 D = 0.03,
                 fib_scale = 1.,
                 det_pos_local = np.array([0., 0., 0.]),
                 head_pos_type = 'r_z_phi',
                 head_position = np.array([1.6, 0. , 90*dtr]),
                 det_dir_type = 'local_theta_phi',
                 det_direction = np.array([50*dtr, 0.*dtr]),
                 tracker = None,
                 arm_rotation = 0.,
                 bundle_fname = 'bundle.npz',
                 comment = '',
                 color = 'r',
                 zero_at_coll = False,
                 ignore_detector_field = False
                 ):
        # this is for the generation of a bundle, calculate a transformation matrix that moves
        # the original z-direction into the direction of the initial velocity
        # detector number e.g. channel
        self.number = number
        # detector name
        self.name = name
        # comment
        self.comment = comment
        # if true set the z-coordinate in the detector coordinate system to 0 after track generation, in that way the detector positions are
        # collimator positions. If False they are detector positions
        self.zero_coll = zero_at_coll
        # Detector Radius
        self.R_det = R_det
        # Collimator Radius
        self.R_coll = R_coll
        # Cylinder radius joining detector to collimator
        self.R_cyl = R_cyl
        # Distance detector-collimator
        self.D = D
        # ignore B-file at detector location
        self.ignore_detector_field = ignore_detector_field
        # scale of angular range for fibbonaci grid (enlarge to make sure with B-fiels all possible trajectories are covered in collimator)
        self.fib_scale = fib_scale
        # detector head position type
        self.head_pos_type = head_pos_type
        # detector head position
        self.head_position = head_position
        # detector direction type
        self.det_dir_type = det_dir_type
        # detector direction 
        self.det_direction = det_direction
        # position and direction setup
        self.det_pos = Cs.pos_dir(pos_type = head_pos_type, 
                                     position = head_position, 
                                     dir_type = det_dir_type, 
                                     direction = det_direction,
                                     position_local = det_pos_local,
                                     alpha = arm_rotation)
        # position in local coord system
        self.det_pos_local = det_pos_local
        # arm rotation
        self.arm_rotation = arm_rotation
        # setup transformation
        self.det_pos.set_detector_to_std()
        # calc position and direction in std coordinate system
        self.det_pos_std = self.det_pos.local_to_std_with_rotation(self.det_pos_local) + self.det_pos.pos
        self.det_dir_std = self.det_pos.local_to_std_with_rotation(self.det_pos.dir_local)
        # set tracker
        if tracker is not None:
            self.set_tracker(tracker)
        self.bundle_fname = bundle_fname
        self.color = color    # color for drawing trajectories
        self.bundle_vars = {'x':0, 'y':1, 'z':2, 'r':3}
        self.Bf_bundle_vars = {'b_pol_r':0, 'b_pol_z':1, 'b_phi':2, 'b_total':3, 'psi_rel':4}

    def set_tracker(self, Tr_l):
        # set the magnetic field tracker used
        self.Tr = Tr_l
        # set initial velocity
        self.v0 = Tr_l.tracker.vmag
        # set particle mass in electron masses
        self.Tr.tracker.particle_mass_me = Tr_l.tracker.particle_mass/Tr_l.constants_and_masses_mod.me
        # magnetic field at detector location
        if self.ignore_detector_field:
            # ignore field
            self.B_det = np.array([0.,0.,0.])
        else:
            # calculate field and use it
            self.B_det = self.Tr.em_fields_mod.bfield3(self.det_pos.pos)

    def init_trajectories(self, N_pos = 10, N_dir = 50, N_s = 50, dN_s = 2):
        # calculatetrajectory initial positions and velocities
        # setup detector bundle get the kinematic parameters from the tracker
        self.bd = MTb.bundle(  N_pos = N_pos,  # number of fibonacci points in detector area
                         N_dir = N_dir, # number of direction fibonacci points
                         N_s = N_s, dN_s = dN_s,        # number of steps inside detector
                         v = self.v0,                    # particle velocity
                         m_me = self.Tr.tracker.particle_mass_me,
                         q_ec = self.Tr.tracker.particle_charge_ec,
                         R_det = self.R_det, R_coll = self.R_coll, R_cyl = self.R_cyl, D = self.D, # detector geometry
                         Bf = self.B_det,                           #  mag. field at detector location
                         scale = self.fib_scale,
                         z_zero = self.zero_coll)
        # intialize
        self.bd.initialize()
        
    def calc_central(self):
        # calculate central trajectory
        # dwetector position
        r_c =  self.det_pos_std          
        # initial velocity along detector-collimator axis
        v_c = self.det_dir_std*self.v0   
        nc = self.Tr.tracker.get_trajectory(r_c, v_c )
        x = C.copy(self.Tr.tracker.trajectory[:nc-1,0])
        y = C.copy(self.Tr.tracker.trajectory[:nc-1,1])
        z = C.copy(self.Tr.tracker.trajectory[:nc-1,2])
        r = np.sqrt(x**2 + y**2)
        self.central_track = np.stack([x,y,z,r]).T
        self.central_bfields = np.array([self.Tr.em_fields_mod.bfield(rr,zz) for rr,zz in zip(r,z)])
        # save additional trajectory bundle  information including kinematics
        general_info = {'particle_charge_ec': self.Tr.tracker.particle_charge_ec+0.,
                      'particle_mass_amu': self.Tr.tracker.particle_mass_amu+0.,
                      'particle_energy_mev': self.Tr.tracker.particle_energy_mev+0.,
                      'efit_file_name':'g'+np.string_(self.Tr.tracker.efit_gfile_name).decode('utf-8').strip(),
                      'r_maxis':self.Tr.flux_par_mod.rmaxis+0.,
                      'z_maxis':self.Tr.flux_par_mod.zmaxis+0.,
                      'particle_velocity':self.v0+0.,
                      'step_size':self.Tr.tracker.step_size+0.,
                      'bundle_variables':self.bundle_vars,
                      'Bf_bundle_variables':self.Bf_bundle_vars,
                      'comment':self.comment }
        np.savez_compressed('central_'+self.bundle_fname,
                            trajectories = self.central_track,
                            B_fields  = self.central_bfields,
                            information = general_info)
        
        
    def calc_trajectories(self, bundle_type = 'full'):
        # ready to track bundles
        t_start = time.time()  # for timing
        # calculate initial values for bundles
        self.bd.calc_bundle(bundle_type = bundle_type)
        # calc trajectory initial positions and acceptance
        r, v, self.acc = self.bd.get_bundle()
        self.r_init_det = r
        self.v_init_det = v
        # transform position to local coordinate system excl arm rotation
        # add local offset due to position, this should move the position to the correct location in the local system
        r_l = (self.det_pos.detector_to_local_no_rotation(r.T)).T + self.det_pos_local
        # transform to Std coordinate system including arm rotation
        r_d =  (self.det_pos.local_to_std_with_rotation(r_l.T)).T
        # transform velocities
        self.v_init = (self.det_pos.detector_to_std(v.T)).T
        # initial positions in tokamak
        self.r_init = r_d + self.det_pos.pos
        # calculate each track and
        # store results in bundle
        bundle = []
        Bf_bundle = [] # magnetic fields along the bundle
        for i,vv in enumerate(self.v_init):
            rr = self.r_init[i]
            # calculate track result is in Tt.tracker.trajectory array
            nc = self.Tr.tracker.get_trajectory(rr, vv )
            # print ('trajectory contains ', nc, ' steps' )
            # need to copy the result values to store them
            x = C.copy(self.Tr.tracker.trajectory[:nc-1,0])
            y = C.copy(self.Tr.tracker.trajectory[:nc-1,1])
            z = C.copy(self.Tr.tracker.trajectory[:nc-1,2])
            r = np.sqrt(x**2 + y**2)
            # calculate the magnetic field used in the calculation
            bf = np.array([self.Tr.em_fields_mod.bfield(rr,zz) for rr,zz in zip(r,z)])
            # combine magnetic field and position
            bundle.append(np.stack([x,y,z,r]).T)
            Bf_bundle.append(bf)
            # all done
        t_end = time.time()
        print("Time used for ", i+1, " tracks = ", t_end - t_start)
        # save trajectory data
        self.bundle = np.array(bundle, dtype = object)
        # save the magetic field
        self.Bf_bundle = np.array(Bf_bundle, dtype = object)
        # save additional trajectory bundle  information including kinematics
        general_info = {'particle_charge_ec': self.Tr.tracker.particle_charge_ec+0.,
                      'particle_mass_amu': self.Tr.tracker.particle_mass_amu+0.,
                      'particle_energy_mev': self.Tr.tracker.particle_energy_mev+0.,
                      'efit_file_name':'g'+np.string_(self.Tr.tracker.efit_gfile_name).decode('utf-8').strip(),
                      'r_maxis':self.Tr.flux_par_mod.rmaxis+0.,
                      'z_maxis':self.Tr.flux_par_mod.zmaxis+0.,
                      'particle_velocity':self.v0+0.,
                      'step_size':self.Tr.tracker.step_size+0.,
                      'bundle_variables':self.bundle_vars,
                      'Bf_bundle_variables':self.Bf_bundle_vars,
                      'comment':self.comment }
        np.savez_compressed(self.bundle_fname,
                            trajectories = self.bundle,
                            B_fields  = self.Bf_bundle,
                            acceptance = self.acc,
                            information = general_info)
