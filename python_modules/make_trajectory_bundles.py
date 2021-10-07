#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 07:05:43 2021

Create trajectory bundles for cylindrical detector-collimator systems

The trajectories initial values are evenly distributed over the detector area
and direction using Fibonacci grids on circle and sphere

Evective acceptance values per trajectory are also calculated

The magnetic field between detector and collimator is assumed to be constant.

The detector-collimator system is assumed to consist of a cylinder with a hole at the beginning
the size of the detector area and a hole at the end the size of the collimator. If both are the same size
as the cylinder the collimator is simply a drilled hole.


|-------------------------------------------|
|                                           |
                         R_cyl              |
R_det                                       R_coll  
                                            
                                            |
|                                           |
|-------------------------------------------|


There are 3 bundle types:
    
There are 3 bundle types:
    
full : - calculate directions and positions on a Fibonacci grid and calculates 
         collimator-detector overlap.
       - includes the effect of the magnetic field in the collimator
       - most realistic but time-consuming calculation

square: - assume square detector and collimator with the same area as a round 
          detector and collimator (R_det and R_coll define the area)
       - ignores the magnetic field in collimator
       - positions and directions on rectangular grids
       - for quick studies

round: - calculates directions on a Fibonacci grid and calculates 
         collimator-detector overlap.
       - start all trajectories at the center of the collimator.
       - ignores the magnetic field in collimator


@author: boeglinw
"""

import numpy as np

# generate fibonacci grids
import fibonacci as FI

# acceptance calculations
import acceptance_array as AA

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
                 R_det = .001, R_coll = .001, D = .05, R_cyl = 0.001,
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
            number of steps for tracking in the magnetic field between detector and collimator. The default is 50

        dN_s: int optional
            additional number of steps for tracking in magnetic field to make sure the entire length of the detector system, is transversed

        v : float , optional
            magnitude of initial velocity. The default is 1e5

        m_me: float, optional
            particle mass in electron masses. The default is 1836.152 (proton mass)

        q_ec: float, optional
            particle charge in elementary charge. The default is 1.

        Bf : np.array(3), float
            constant magnetic field between deteactor and collimator. The default is np.array([0., 0., 0.]

        R_det : float , optional
            detector radius. The default is .001.

        R_coll : float, optional
            collimator radis. The default is .001.
            
        R_cyl : float, optional
            cylinder radius joining detector and collimaotr 

        D : float, optional
            distance detector-collimator. The default is .05.

        Ef : np.array(3), float
            electric field, The default is np.array([0., 0., 0.]

        scale: float, optional
            scale factor to increase maximal polar direction angle. The default is 1.

        z_zero: bool, optional
            set z-coordinate of final position to 0. This is useful if detector positions are referring to collimator position. The default is False

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
        self.R_cyl = R_cyl
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
        
    def calc_bundle(self, stop = True, bundle_type = 'full'):
        # calculate selected bundle type
        if bundle_type == 'full':
            self.calc_full_bundle(stop = stop)
        elif bundle_type == 'square':
            self.calc_square_bundle()
        elif bundle_type == 'round':
            self.calc_round_bundle()
        else:
            print(f'Unknown bundle type : {bundle_type}, nothing calculated !')
        return

    def calc_full_bundle(self, stop = True):
        # initialize arrays
        results = []
        t_det = []
        t_coll = []
        hits = []
        # calculate trajectories
        #  stop when wall is hit
        BC.boris_cylinder.stop_at_hits = stop
        # calculate trajectory for each position and each direction
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
                BC.boris_cylinder.init(self.q_ec, self.m_me, self.R_cyl, self.R_coll, self.D, step_l, Nstep_l, vs )
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


    def calc_square_bundle(self):
        """
        Make a trajectory bundle with square collimator and detector. The area of the square is the same
        as the area of a circle with the given radius
    
        Returns
        -------
        None.
    
        """
        s_det = np.sqrt(np.pi*self.R_det**2)/2.
        s_col = np.sqrt(np.pi*self.R_coll**2)/2.
        # make detector grid
        # N_pos = numner of detecor grid points
        # N_dir = number of collimator grid points
        dx_d = s_det/self.N_pos
        dx_c = s_col/self.N_dir
        # 1d grid positions 
        x_d = (np.arange(self.N_pos) + 0.5)*dx_d - s_det
        y_d = x_d
        x_c = (np.arange(self.N_dir) + 0.5)*dx_c - s_col
        y_c = x_c
        # detector grid positions
        xxd, yyd = np.meshgrid(x_d, y_d)
        # collimator grid positions
        xxc, yyc = np.meshgrid(x_c, y_c)
        # detector collimator combinations
        XD,XC = np.meshgrid(xxd.flatten(), xxc.flatten())
        xd = XD.flatten(); xc = XC.flatten()
        YD,YC = np.meshgrid(yyd.flatten(), yyc.flatten())
        yd = YD.flatten(); yc = YC.flatten()
        # calculate trajectories
        zc = np.ones_like(xd)*self.D
        vx = xc - xd
        vy = yc - yd
        vz = zc
        #
        V = np.array([vx,vy,vz])
        v_mag = np.apply_along_axis(np.linalg.norm, 0, V)
        # velocity unit vectors
        uv = V/v_mag
        # initial positions
        self.r_coll = (np.array([xc, yc, np.zeros_like(xc)])).T
        # initical velocity
        self.v_coll = (uv*self.v).T
        # calc acceptance
        acc_x = AA.accept1d(s_col/self.N_dir, s_det/self.N_pos, self.D, (xd - xc))
        acc_y = AA.accept1d(s_col/self.N_dir, s_det/self.N_pos, self.D, (yd - yc))
        self.dAcc = acc_x*acc_y
        self.ratio = 1.
        

    def calc_round_bundle(self):
        """
        Make a trajectory bundle with round collimator and detector. Using a fibonacci grid for the angles.
    
        Returns
        -------
        None.
    
        """
        fg = FI.fibonacci_grid(self.N_dir)
        
        # largest possible angle
        th_max = np.arctan((self.R_det + self.R_coll)/self.D)
        
        # total solid angle
        A = 2.*np.pi*(1. - np.cos(th_max))
        # solid angle per point
        dA = A/self.N_dir
        
        phi, theta = fg.get_sphere(0., th_max)
        # initial directions
        xp = np.cos(phi)*np.sin(theta)
        yp = np.sin(phi)*np.sin(theta)
        zp = np.cos(theta)
        # velocity unit vector
        self.v_coll = (np.array([xp, yp, zp])*self.v).T
        # collimator positions
        xc = np.zeros_like(xp)
        self.r_coll = (np.array([xc, xc, xc])).T
        # calculate overlapp
        over = AA.C_overlap(theta, self.R_coll, self.R_det, self.D)
        # acceptance per trajectory
        self.dAcc = over*np.cos(theta)*dA
        self.ratio = 1.

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


