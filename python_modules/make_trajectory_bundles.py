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
                 v = 1.38e7,
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
            magnitude of initial velocity. The default is 1.38e7 (~ 1 MeV)

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
        self.A_coll = np.pi * self.R_coll**2
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
        """
        Calculate a set of initial values for a number of direction and collimator position points 

        Parameters
        ----------
        stop : TYPE, optional
            DESCRIPTION. The default is True.
            stop the calculation of the trajector between detector and collimator if the cylinder wall is hit

        Returns
        -------
        None.


        The following quantities are calculated:
            
            self.r_coll 
                initial positions
                
            self.v_coll
                initial velocities
                
            self.dAcc
                array of acceptances for each trajectory


        """
        # initialize arrays
        results = []
        t_det = []
        t_coll = []
        wall_hits = []
        coll_hits = []
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
                # save results track is an array of 6-dim points 
                # tr[:,0:2] : positions
                # tr[:,3:6] : velocities
                tr = CO.copy(BC.boris_cylinder.track[:,:])
                # get initial and final values
                # select only parts inside the system
                sel = tr[:,2]<= self.D
                t_init = tr[sel][0]
                t_final = tr[sel][-1]
                if self.save_track:
                    results.append(tr)  # save position and  velocities
                t_det.append(t_init)
                t_coll.append(t_final)
                # check for a wall hit
                wall_hits.append(CO.copy(BC.boris_cylinder.hit))
                # check for collimator hit
                r_xy_coll = np.sqrt(t_final[0]**2 + t_final[1]**2)
                coll_hits.append(r_xy_coll >= self.R_coll )
                

        print(f'time = {timeit.default_timer() - start_time}')

        # make arrays

        if self.save_track:
            self.results = np.array(results, dtype = object)            
        self.t_det = np.array(t_det)
        self.t_coll = np.array(t_coll)

        self.wall_hits = np.array(wall_hits )> 0
        self.coll_hits = np.array(coll_hits)
        # either wall or coll has been hit
        hits = self.coll_hits | self.wall_hits


        self.ok = ~hits

        # fraction of tracks that hit wall
        self.ratio = np.count_nonzero(self.ok)/hits.shape[0]

        print(f'passing fraction = {self.ratio:.3f}')

        # calculate acceptance for the passed particles
        if self.save_track:
            results_p = self.results[self.ok]

            # particle positions at collimator
            self.r_coll = results_p[:,-1,:3]
            # particle velocities at collimator
            self.v_coll = results_p[:,-1,3:]
            # particle velocity at detector
            self.v_init = results_p[:,0,3:]
        else:
            # particle positions at collimator
            self.r_coll = self.t_coll[self.ok,:3]
            # particle velocities at collimator
            self.v_coll = self.t_coll[self.ok,3:]
            # particle velocity at detector
            self.v_init = self.t_coll[self.ok,3:]
            
        self.v_coll_mean = np.apply_along_axis(np.mean, 0, self.v_coll)


        # calc effective acceptance
        # cos(theta) of final direction
        cth = (self.v_coll/self.v)[:,2]
        cthi = (self.v_init/self.v)[:,2]


        # calculate acceptances
        # acceptance per track
        self.dAcc = cth*self.dA_pos*self.dA_dir
        self.dAcci = cthi*self.dA_pos*self.dA_dir
        # using final values
        self.Acc = np.sum(self.dAcc)
        # using initial values
        self.Acci = np.sum(self.dAcci)

    def calc_round_transmission(self, theta, phi, return_details = False):
        """
        calculate transmission factor for a selected direction for a general detector including 
        magnetic field at detector location. Select a position at the collimator and back-track it
        to the detector to see if it hit the detector

        Returns
        -------
        None.

        """
        # initialize result arrays
        t_results = []
        t_r_det = []
        t_ok = []
        t_ratio = []
        # create collimator position array
        r_coll_init, th_pos = fg.get_circle(self.R_coll, N = self.N_pos) 
        xc = r_coll_init*np.cos(th_pos)
        yc = r_coll_init*np.sin(th_pos)
        # displacement vector at collimator
        delta_r_coll = np.array([xc, yc, np.zeros_like(xc)]).T
        #
        start_time = timeit.default_timer()
        
        # reverse field for back tracking
        BC.boris_cylinder.set_b0(-self.B)
        # do not stop if trajectory hits a wall as it is only needed for the 
        BC.boris_cylinder.stop_at_hits = False  
        Nstep_l = self.N_s + self.dN_s
        # loop over initial velocities
        for k,th in enumerate(theta):
            ph = phi[k]
            results = []
            r_det = []
            r_coll = []
            hits = []
            # calculate trajectory for each position and a given direction
            # inital velocity direction
            vs = -np.array([np.sin(th)*np.cos(ph),
                           np.sin(th)*np.sin(ph),
                           np.cos(th)])*self.v
    
            step_l = self.step/np.cos(th)                        
            # calculate central trajectory  
            # init boris_cylinder for this velocity
            BC.boris_cylinder.init(self.q_ec, self.m_me, self.R_cyl, self.R_det, self.D, step_l, Nstep_l, vs )
            # calc. reverse orbit
            nc0 = BC.boris_cylinder.track_cylinder(np.array([0.,0.,self.D]))
            # store central orbit
            tr0 = CO.copy(BC.boris_cylinder.track[:,:])
            r0 = tr0[:,0:3]  # central trajectory (use positions only)
            # parallel displace central trajectory     
            # loop over all initial collimator positions
            for rs in delta_r_coll:
                # move central trajectory to rs
                rr = rs + r0
                xr = rr[:,0]; yr = rr[:,1];zr = rr[:,2]
                # get initial and final values
                Rr = np.sqrt(xr**2 + yr**2)  # transverse position
                sel = zr >= 0.
                # determine if the cylinder wall was hit or the detector was missed
                hit_w = (Rr[sel] >= self.R_cyl).max()
                miss_det = (Rr[sel] >= self.R_det).max()
                hit = hit_w or miss_det  # record the result
                r_init = rr[sel][0]
                r_final = rr[sel][-1]
                r_det.append(r_final)
                r_coll.append(r_init)
                hits.append(hit)        # save if it hit somewhere
            hits = np.array(hits)
            ok = ~hits            
            # calculate transmission factor
            ratio = np.count_nonzero(ok)/hits.shape[0]
            t_ratio.append(ratio)
        
            if return_details:
                 t_r_det.append(np.array(r_det)) # array of track positions in detector
                 t_ok.append(np.array(ok))  # array of which track position lead to a hit   
        if return_details:
            return np.array(t_ratio), np.array(r_coll), np.array(t_r_det),  np.array(t_ok)
        else:
            return np.array(t_ratio), np.array(r_coll)  # ratio and array of collimator positions
        
        
        
        
        

    def calc_square_bundle(self):
        """
        Make a trajectory bundle with square collimator and detector. The area of the square is the same
        as the area of a circle with the given radius
    
        Returns
        -------
        None.
    
 
        The following quantities are calculated:
            
            self.r_coll 
                initial positions
                
            self.v_coll
                initial velocities
                
            self.dAcc
                array of acceptances for each trajectory


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
        

    def calc_round_bundle(self, full_bundle = False):
        """
        Make a trajectory bundle with round collimator and detector. Using a fibonacci grid for the angles.
        
        Parameters
        ----------
        
        full_bundle: bool (optional) 
        
            True : similar to the output of calc_full_bundle, calcultaes a set of initial values for 
                   trajectory bundles with different positions in the collimator.
                   
            False(default): calculates a set of directions and the associated transmission factors. All initial positions
                            are in the center of the collimator. This allows for an efficient calculation of
                            rates
                      
              
    
        Returns
        -------
        None.
        
        
        The following quantities are calculated:
            
            self.r_coll 
                initial positions
                
            self.v_coll
                initial velocities
                
            self.dAcc
                array of acceptances for each trajectory
        
        
    
        """
        fg = FI.fibonacci_grid(self.N_dir)
        
        # largest possible angle
        th_max = np.arctan((self.R_det + self.R_coll)/self.D)
        
        # total solid angle
        A = 2.*np.pi*(1. - np.cos(th_max))
        # solid angle per point
        dA = A/self.N_dir
        
        # collimator area element
        dA_coll = self.A_coll/self.N_pos
        
        
        phi, theta = fg.get_sphere(0., th_max)
        self.theta_dir = theta
        self.phi_dir = phi
        # initial directions
        xp = np.cos(phi)*np.sin(theta)
        yp = np.sin(phi)*np.sin(theta)
        zp = np.cos(theta)
        # velocity unit vectors
        v_coll = (np.array([xp, yp, zp])*self.v).T

        # calculate overlapp for round bundle and check if it is ok incl B-field
        if full_bundle :
            ratio, r_coll, r_det, ok_a = self.calc_round_transmission(theta, phi, return_details = True)
            self.c_ok = ok_a
        else:
            ratio, r_coll = self.calc_round_transmission(theta, phi, return_details = False)
        self.c_ratio = ratio   # for debugging
        self.c_r_coll = r_coll
        
        # if a full bundle is needed store the collimator positions (this should be identical to the full bundle)
        if full_bundle:
            self.r_coll = np.vstack([r_coll[ok] for ok in ok_a])
            # store the 
            self.v_init_l = []
            for i,th in enumerate(theta):
                n_ok = r_coll[ok_a[i]].shape[0]
                if n_ok == 0:
                    self.v_init_l.append( np.empty(shape = (0,3)) )
                else:
                    self.v_init_l.append(n_ok * [self.v_coll[0]])
            self.v_coll = np.vstack(self.v_init_l)
            self.dAcc = np.concatenate( [dA_coll*np.cos(th)*dA*np.ones_like(r_coll[ok_a[i]][:,0]) for i,th in enumerate(theta)])
        else:
            # set the positions at the center of the collimator
            xc = np.zeros_like(xp)
            self.r_coll = (np.array([xc, xc, xc])).T
            self.dAcc = ratio*self.A_coll*np.cos(theta)*dA
            self.ratio = ratio
            self.v_coll = v_coll
        # acceptance per trajectory

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


