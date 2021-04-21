#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 07:05:43 2021

Create trajectory bundles for cylindrical detector-collimator systems

The trajectories initial values are evenly distributed over the detector area
and direction using Fibonacci grids on circle and sphere

Evective acceptance values per trajectory are also calculated

@author: boeglinw
"""

import numpy as np

# generate fibonacci grids
import fibonacci as FI

# calc culate trjectories in constant fields
import BorisCylpy as BC
import copy as CO

# timing purposes
import timeit

dtr = np.pi/180.

# initialize
fg = FI.fibonacci_grid()


class bundle:
    
    def __init__(self,  N_pos = 50, N_dir = 100, 
                 N_s = 50, dN_s = 5, 
                 v = 1.e5,
                 m_me = 1836.152,
                 q_ec = 1., 
                 Bf = np.array([0., 0., 0.]), 
                 R_det = .001, R_coll = .001, D = .05, 
                 Ef = np.array([0., 0., 0.]), 
                 scale = 1., 
                 z_zero = False):
        """
        
        calculate the initial values for a bundle of trajectories evenly distributed over the 
        acceptance of a circular detector-collimator system.

        Parameters
        ----------
        N_pos : int , optional
            number of fibonacci points for initial positions The default is 50
            
        N_dir: int , optional
            number of fibonacci points for directions The default is 100
            
        N_s: int optional
            number of steps for tracking in magn. field between detector and collimator. The default is 50

        dN_s: int optional
            additional number of steps for tracking in magn. field to make sure the entire length of the detector system, is transversed
            
        v : float , optional
            magnitude of initial velocity. The default is 1e5

        m_me: float, optional
            particle mass in electron masses. The default is 1836.152 (proton mass)
            
        q_ec: float, optional
            particle charge in elementary charge. The default is 1.
            
        Bf : np.array(3), float
            magnetic field. The default is np.array([0., 0., 0.]
                                                    
        R_det : float , optional
            detector radius. The default is .001.
            
        R_coll : float, optional
            collimator radis. The default is .001.
            
        D : float, optional
            distance detector-collimator. The default is .05.
            
        Ef : np.array(3), float
            electric field, The default is np.array([0., 0., 0.]
                                                    
        scale: float, optional
            scale factor to increase maximal polar direction angle
            
        z_zero: bool, optional
            set z-coordinate of final position to 0. This is usefule of detector positions are referring to collimator position.
            
        Returns
        -------
        None.

        """
        self.N_pos = N_pos
        self.N_dir = N_dir
        self.N_s = N_s
        self.v = v
        self.q_ec = q_ec
        self.m_me = m_me
        self.B = Bf
        self.Ef = Ef
        self.R_det = R_det
        self.R_coll = R_coll
        self.D = D
        self.scale = scale
        # additional number of boris steps to make sure the entire length of the detector 
        # has been traversed
        self.dN_s = dN_s
        self.save_track = False
        self.z_zero = z_zero
        
        
        
    def initialize(self):
        """
        Initialize the calculation, this should always be done if one of the parameters has changed

        Returns
        -------
        None.

        """
        # setup mag. field tracing
        BC.boris_cylinder.set_b0(self.B)
        self.step = self.D/self.N_s
        # initialize calc
            
        # sunflower r theta distribution for initial detector position
        self.r_init, self.th_pos = fg.get_circle(self.R_det, N = self.N_pos)
        self.A_det = np.pi * self.R_det**2
        self.dA_pos = self.A_det/self.N_pos
        
        # sunflower for directions
        # largest angle
        self.th_max = np.arctan((self.R_det + self.R_coll)/self.D) *self.scale
        # total solid angle
        self.A_dir = 2.*np.pi*(1. - np.cos(self.th_max))
        # solid angle per point
        self.dA_dir = self.A_dir/self.N_dir
        self.phi_dir, self.theta_dir =  fg.get_sphere(0., self.th_max, N = self.N_dir)
        
    def calc_bundle(self, stop = True):
        # initialize arrays
        results = []
        t_det = []
        t_coll = []
        hits = []
        # calculate trajectories
        #  stop when wall is hit
        BC.boris_cylinder.stop_at_hits = stop
        # calculate trajectory for reach position and each direction
        start_time = timeit.default_timer()
        for i,rr in enumerate(self.r_init):
            # initial position
            x = rr*np.cos(self.th_pos[i])
            y = rr*np.sin(self.th_pos[i])
            rs = np.array([x, y, 0.])
            for j, phv in enumerate(self.phi_dir):
                thv = self.theta_dir[j]
                # inital velocity
                vs = np.array([np.sin(thv)*np.cos(phv),
                               np.sin(thv)*np.sin(phv),
                               np.cos(thv)])*self.v
                # initialize calc
                # number of boris steps
                Nstep_l = self.N_s + self.dN_s
                # make sure step size and Nstep allow for a trajectory to travel the full length (d) 
                # of the detector
                # 
                step_l = self.step/np.cos(thv)
                BC.boris_cylinder.init(self.q_ec, self.m_me, self.R_det, self.R_coll, self.D, step_l, Nstep_l, vs )
                # calc trajectory
                nc = BC.boris_cylinder.track_cylinder(rs)
                # save results
                tr = CO.copy(BC.boris_cylinder.track[:,:])
                # get initial and final values
                # select only parts inside the system

                sel = tr[:,2]<= self.D
                t_init = tr[sel][0]
                t_final = tr[sel][-1]
                                 
                if self.save_track:
                    results.append(tr)  # save position and  velocities
                else:
                    t_det.append(t_init)
                    t_coll.append(t_final)
                hits.append(CO.copy(BC.boris_cylinder.hit))
        
        print(f'time = {timeit.default_timer() - start_time}')
        
        
        # make arrays
        
        hits = np.array(hits)>0
        if self.save_track:
            self.results = np.array(results, dtype = object)
        else:
            self.t_det = np.array(t_det)
            self.t_coll = np.array(t_coll)
        
        self.ok = ~hits            
        
        # fraction of tracks that hit wall
        self.ratio = np.count_nonzero(self.ok)/hits.shape[0]
        
        print(f'passing fraction = {self.ratio:.3f}')
        
        # calculate acceptance for the passed paricles
        if self.save_track:
            results_p = self.results[self.ok]
            
            # particle positions at collimator
            self.r_coll = results_p[:,-1,:3]            
            # particle velocities at collimator
            self.v_coll = results_p[:,-1,3:]
        else:
            # particle positions at collimator
            self.r_coll = self.t_coll[self.ok,:3]            
            # particle velocities at collimator
            self.v_coll = self.t_coll[self.ok,3:]
            
        self.v_coll_mean = np.apply_along_axis(np.mean, 0, self.v_coll)
        
        
        # calc effective acceptance
        # cos(theta) of final direction
        cth = (self.v_coll/self.v)[:,2]
        cthi = (self.v_coll/self.v)[:,2]
         
        
        # calculate acceptances
        # acceptance per track
        self.dAcc = cth*self.dA_pos*self.dA_dir
        self.dAcci = cthi*self.dA_pos*self.dA_dir
        # using final values 
        self.Acc = np.sum(self.dAcc)
        # using initial values
        self.Acci = np.sum(self.dAcci)


    def get_bundle(self):
        """
        Return position, velocities and acceptances for a bundle of trajectories. These can be used later in
        an orbit program.

        Returns
        -------
        r : np.array(:,3), float
            initial positions.
            
        v : np.array(:,3), float
            initial velocities.
            
        A : np.array(:), float
            detector acceptance for each track.

        """
        r,v, A  = self.r_coll, self.v_coll, self.dAcc
        
        # force z-component of r to be not larger than D
        
        if self.z_zero:
            # set z-coordinate of initial position to 0, i.e. coordinate system origin is collimator
            r[:,2]-= r[:,2]
        return r,v,A 
    
    
