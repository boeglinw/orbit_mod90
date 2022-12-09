#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:56:52 2021

function to load a trajectory bundle

usage example::

    >>> import numpy as np
    >>> import matplotlib.pyplot as pl
    >>> import load_bundle as LB

    >>> b0 = LB.trajectory_bundle('../example_data/det_0_Channel0.npz')

    # loop over all trajectories and make ab r-z plot

    >>> for tr in b0.trajectories():
            x,y,z,r = tr
            pl.plot(r,z)
    >>> pl.gca().set_aspect('equal')

@author: boeglinw
"""


import numpy as np


class trajectory_bundle:

    def __init__(self, fname):
        """
        Load trajectory bundle

        Parameters
        ----------
        fname : String
            file name.

        Returns
        -------
        trajectory bundle object.

        """
        self.filename = fname
        self.load_data()

    def load_data(self):
        """
        Load the data stored in self.filename

        Returns
        -------
        None.

        """
        d = np.load(self.filename, allow_pickle = True)
        self.keys = list(d.keys())
        self.bundle = d['trajectories']
        self.n_trajectories = len(self.bundle)
        self.Bf = d['B_fields']
        self.acceptance = d['acceptance']
        self.information = d['information'].item()

    def trajectories(self):
        """
        Generator to access all trajectpries

        Yields
        ------
        generator
            to loop over stored trajectories.

        # loop over all trajectories and make ab r-z plot

        >>> for tr in b0.trajectories():
                x,y,z,r = tr

        """
        # generator for the trajectories in the bundle
        counter = 0
        while counter < self.n_trajectories:
            yield np.array(list((self.bundle[counter]).T))
            counter += 1

    def B_fields(self):
        """
        Generator to access the B-fields along the trajectories

        Yields
        ------
        generator
            to loop over the B-field along each trajectory.

        """
        # generator for the B-fields along the trajectory
        counter = 0
        while counter < self.n_trajectories:
            yield (self.Bf[counter]).T
            counter += 1


    def get_trajectory(self, i, component = None):
        """
        Returns trajectory number i as a numpy array

        Parameters
        ----------
        i : int
            trajectory number.
        component: string, optional
            name of position vector component: x, y, z, r
            

        Returns
        -------
        np.array
            trajectory as np.array([x,y,z,r]).

        """
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        
        if component is None:
            return np.array(list((self.bundle[i]).T))
        elif component in self.information['bundle_variables'].keys():
            cc = self.information['bundle_variables'][component]
            return (np.array(list((self.bundle[i]).T)))[cc,:]
        else:
            print(f"No such component : {component}, possible values {list(self.information['bundle_variables'].keys())}")
            return None        


    def get_B_fields(self, i, component = None, cartesian = False):
        """
        Returns the magnetic fields along trajectory number i as a numpy array

        Parameters
        ----------
        i : int
            trajectory number.
        
        component: string, optional
            name of position vector component: b_pol_r, b_pol_z, b_phi, b_total, psi_rel

        Returns
        -------
        np.array
            B-field along trajectory.

        """
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        if component is None:
            if cartesian:
                # calc unit vectors
                x,y,z,r = list(self.get_trajectory(i))
                er = np.array([x/r, y/r, np.zeros_like(x)])
                ephi = np.array([-er[1], er[0], np.zeros_like(x)] )
                ez = np.array([np.zeros_like(x), np.zeros_like(x), np.ones_like(x)])
                # get field comnponents
                br, bz, bphi, bt, psi = list(self.Bf[i].T)
                # calc cartesian field vecgor
                bf = br*er + bz*ez + bphi*ephi
                return bf
            else:
                return self.Bf[i].T
        elif component in self.information['Bf_bundle_variables'].keys():
            cc = self.information['Bf_bundle_variables'][component]
            return (self.Bf[i].T)[cc,:]
        else:
            print(f"No such component : {component}, possible values {list(self.information['Bf_bundle_variables'].keys())}")
            return None
    
    def __len__(self):
        """
        return number of trajectories
        """
        return self.n_trajectories

    def __getitem__(self,i):
        """
        return trajectory number i
        """
        return (self.bundle[i]).T


