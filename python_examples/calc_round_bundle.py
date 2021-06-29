#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 07:11:53 2021

@author: boeglinw
"""
import numpy as np

import acceptance_array as AA
import fibonacci as FI


def calc_round_bundle(R_det, R_col, N_dir, D):
    """
    Make a trajectory bundle with round collimator and detector. Using a fibonacci grid for the angles.

    Returns
    -------
    None.

    """
    fg = FI.fibonacci_grid(N_dir)
    
    #indices = np.arange(0, num_pts, dtype=float) + 0.5
    #arg = 1 - 2*indices/num_pts
    
    # largest angle
    th_max = np.arctan((R_det + R_col)/D)
    
    # total solid angle
    A = 2.*np.pi*(1. - np.cos(th_max))
    # solid angle per point
    dA = A/N_dir
    
    phi, theta = fg.get_sphere(0., th_max)
    # initial directions
    xp = np.cos(phi)*np.sin(theta)
    yp = np.sin(phi)*np.sin(theta)
    zp = np.cos(theta)
    # velocity unit vector
    uv = np.array([xp, yp, zp])
    # collimator positions
    xc = np.zeros_like(xp)
    r_col = np.array([xc, xc, xc])
    # calculate overlapp
    over = AA.C_overlap(theta, R_col, R_det, D)
    # acceptance per trajectory
    accept = over*np.cos(theta)*dA
    return r_col, uv, accept