#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 07:11:53 2021

@author: boeglinw
"""
import numpy as np

import acceptance_array as AA


def calc_square_bundle(R_det, N_pos, R_col, N_dir, D):
    """
    Make a trajectory bundle with square collimator and detector. The area of the square is the same
    as the area of a circle with the given radius

    Returns
    -------
    None.

    """
    s_det = np.sqrt(np.pi*R_det**2)/2.
    s_col = np.sqrt(np.pi*R_col**2)/2.
    # make detector grid
    # N_pos = numner of detecor grid points
    # N_dir = number of collimator grid points
    dx_d = s_det/N_pos
    dx_c = s_col/N_dir
    # 1d grid positions 
    x_d = (np.arange(N_pos) + 0.5)*dx_d - s_det
    y_d = x_d
    x_c = (np.arange(N_dir) + 0.5)*dx_c - s_col
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
    zc = np.ones_like(xd)*D
    vx = xc - xd
    vy = yc - yd
    vz = zc
    #
    V = np.array([vx,vy,vz])
    v_mag = np.apply_along_axis(np.linalg.norm, 0, V)
    # velocity unit vectors
    uv = V/v_mag
    # initial positions
    r_col = np.array([xc, yc, zc])
    # calc acceptance
    acc_x = AA.accept1d(s_col/N_dir, s_det/N_pos, D, (xd - xc))
    acc_y = AA.accept1d(s_col/N_dir, s_det/N_pos, D, (yd - yc))
    accept = acc_x*acc_y
    return r_col, uv, accept