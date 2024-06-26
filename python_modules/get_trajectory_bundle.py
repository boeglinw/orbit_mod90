#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:56:52 2021

function to load a trajectory bundle

usage example::

    >>> import numpy as np
    >>> import matplotlib.pyplot as pl
    >>> import get_trajectory_bundle as GTB

    >>> b0 = GTB.trajectory_bundle('../example_data/det_0_Channel0.npz')

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
            file name containing the bundle. This should have been written by detectors.calc_trajectories

        Returns
        -------
        trajectory bundle object.

        >>> b0 = GTB.trajectory_bundle('../example_data/det_0_Channel0.npz')

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
        self.kinematics = d['kinematics'].item()

    def trajectories(self):
        """
        Generator to access all trajectpries

        Yields
        ------
        generator
            to loop over stored trajectories.

        loop over all trajectories and make ab r-z plot:

        >>> for tr in b0.trajectories():
                x,y,z,r = tr
                pl.plot(r,z)
        >>> pl.gca().set_aspect('equal')

        """
        # generator for the trajectories in the bundle
        counter = 0
        while counter < self.n_trajectories:
            yield np.array(list(self.bundle[counter]))
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
            yield self.Bf[counter]
            counter += 1


    def get_trajectory(self, i):
        """
        Returns trajectory number i as a numpy array

        Parameters
        ----------
        i : int
            index.

        Returns
        -------
        np.array
            trajectory number i as np.array([x,y,z,r]).

        """
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        return np.array(list(self.bundle[i]))

    def get_B_fields(self, i):
        """
        Returns the magnetic fields along trajectory number i as a numpy array

        Parameters
        ----------
        i : int
            trajectory number

        Returns
        -------
        np.array
            B-field along trajectory number i

        """
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        return self.Bf[i].T

    def __len__(self):
        """
        return number of trajectories:

        >>> len(b0)   # returns the number of trajectories in bundle b0

        """
        return self.n_trajectories

    def __getitem__(self,i):
        """
        return trajectory number i

        >>> b0[i] # returns trajectory number i

        """
        return self.bundle[i]


