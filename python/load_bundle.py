#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:56:52 2021

function to load a trajectory bundle

usage example:
    
import numpy as np
import matplotlib.pyplot as pl
import load_bundle as LB
 
b0 = LB.trajectory_bundle('../example_data/det_0_Channel0.npz')

# loop over all trajectories and make ab r-z plot

for tr in b0.trajectories():
    x,y,z,r = tr
    pl.plot(r,z)
pl.gca().set_aspect('equal')

@author: boeglinw
"""


import numpy as np


class trajectory_bundle:
    def __init__(self, fname):
        self.filename = fname
        self.load_data()
        
    def load_data(self):
        d = np.load(self.filename, allow_pickle = True)
        self.keys = list(d.keys())
        self.bundle = d['trajectories']
        self.n_trajectories = len(self.bundle)
        self.Bf = d['B_fields']
        self.acceptance = d['acceptance']
        self.kinematics = d['kinematics'].item()
        
        
    def trajectories(self):
        # generator for the trajectories in the bundle
        counter = 0
        while counter < self.n_trajectories:
            yield np.array(list(self.bundle[counter]))
            counter += 1
            
    def B_fields(self):
        # generator for the B-fields along the trajectory
        counter = 0
        while counter < self.n_trajectories:
            yield self.Bf[counter]
            counter += 1
            

    def get_trajectory(self, i):
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        return np.array(list(self.bundle[i]))

    def get_B_fields(self, i):
        if i >= self.n_trajectories :
            print( f'I have only {self.n_trajectories} trajectories, number {i} does not exist')
            return None
        return self.Bf[i].T
