#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 07:11:37 2021

@author: boeglinw
"""
import copy as C
import time
import numpy as np


# coordinate transformations
import coordinate_systems as Cs
# detector-collimator trajectory setup
import trajectory_bundles as Tb


dtr = np.pi/180.


#%% detector class

class detector:
    def __init__(self, number,  name,
                 R_det = 0.001,
                 R_coll = 0.001,
                 D = 0.03,
                 pos_local = np.array([0., 0., 0.]),
                 pos_type = 'r_z_phi',
                 position = np.array([1.6, 0. , 90*dtr]),
                 dir_type = 'local_theta_phi',
                 direction = np.array([50*dtr, 0.*dtr]),
                 tracker = None,
                 rotation = 0.,
                 bundle_fname = 'bundle.npz',
                 comment = ''
                 ):
        # this is for the generation of a bundle, calculate a trasnformation matrix that moves
        # the original z-direction into the direction of the initial velocity
        self.number = number
        self.name = name
        self.comment = comment
        self.R_det = R_det
        self.R_coll = R_coll
        self.D = D
        # position and direction setup
        self.det_pos = Cs.pos_dir(pos_type = pos_type, position = position, dir_type = dir_type, direction = direction, alpha = rotation)   
        # position in local coord system
        self.pos_local = pos_local
        # arm rotation
        self.rotation = rotation
        # setup transformation
        self.det_pos.set_detector_to_std()
        if tracker is not None:
            self.set_tracker(tracker)
        self.bundle_fname = bundle_fname
        
    def set_tracker(self, Tr_l):
        # set the magnetic field tracker used
        self.Tr = Tr_l
        # set initial velocity
        self.v0 = Tr_l.tracker.vmag
        # set particle mass in electron masses
        self.Tr.tracker.particle_mass_me = Tr_l.tracker.particle_mass/Tr_l.constants_and_masses_mod.me
        # magnetic field at detector location
        self.B_det = self.Tr.em_fields_mod.bfield3(self.det_pos.pos)           

    def init_trajectories(self, N_pos = 10, N_dir = 50, N_s = 50, dN_s = 2, z_zero = False):                
        # calculatetrajectory initial positions and velocities                     
        # setup detector bundle get the kinematic parameters from the tracker
        self.bd = Tb.bundle(  N_pos = N_pos,  # number of fibonacci points in detector area
                         N_dir = N_dir, # number of direction fibonacci points                                
                         N_s = N_s, dN_s = dN_s,        # number of steps inside detector 
                         v = self.v0,                    # particle velocity
                         m_me = self.Tr.tracker.particle_mass_me,
                         q_ec = self.Tr.tracker.particle_charge_ec,
                         R_det = self.R_det, R_coll = self.R_coll, D = self.D, # detector geometry 
                         Bf = self.B_det,                           #  mag. field at detector location
                         z_zero = z_zero)
        # intialize 
        self.bd.initialize()
        
    def calc_trajectories(self):
        # ready to track bundles
        t_start = time.time()  # for timing
        # calculate initial values for bundles
        self.bd.calc_bundle()        
        # calc trajectory initial positions and acceptance
        r, v, self.acc = self.bd.get_bundle()
        # transform position to local coordinate system excl arm rotation
        # add local offset due to position, this should move the position to the correct location in the local system
        r_l = (self.det_pos.detector_to_local_no_rotation(r.T)).T + self.pos_local
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
            bundle.append((x, y, z, r))  # store trajectory information in a bundle 
            # calculate the magnetic field used in the calculation
            r_pos = np.array([x,y,z]).T
            Bf_bundle.append( np.array([self.Tr.em_fields_mod.bfield3(r_pos_l) for r_pos_l in r_pos]) )
            # all done
        t_end = time.time()
        print("Time used for ", i, " tracks = ", t_end - t_start)
        # save trajectory data
        self.bundle = np.array(bundle, dtype = object)
        # save the magetic field
        self.Bf_bundle = np.array(Bf_bundle, dtype = object)
        # save additional trajectory bundle  information including kinematics
        particle_kin = {'particle_charge_ec': self.Tr.tracker.particle_charge_ec, 
                      'particle_mass_amu': self.Tr.tracker.particle_mass_amu,
                      'particle_energy_mev': self.Tr.tracker.particle_energy_mev,
                      'particle_velocity':self.v0,
                      'comment':self.comment }
        np.savez_compressed(self.bundle_fname, 
                            trajectories = self.bundle, 
                            B_fields  = self.Bf_bundle,
                            acceptance = self.acc, 
                            kinematics = particle_kin)

class trajectory_bundle:
    """
    Class to load a bundle of trajectories beloning to a detector
    """
    def __init__(self, filename):
        """
        create an instance and load the file

        Parameters
        ----------
        filename : string
            trajectory bundle filename (npz file_.

        Returns
        -------
        None.

        """
        self.filename = filename
        self.load_trajectories()
            
    def load_trajectories(self):
        """
        load the npz file and store the data

        Returns
        -------
        None.

        """
        with np.load(self.filename, allow_pickle=True) as t_data:
            self.trajectories = t_data['trajectories']
            self.acceptance = t_data['acceptance']
            self.B_fields =t_data['B_fields']
            # get the kinematics dictionary
            self.kinematics = t_data['kinematics'].item()
            self.number_of_trajectories = len(self.trajectories)
            

            
            
                
                