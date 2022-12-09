#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 10:13:12 2022

@author: boeglinw
"""

import numpy as np
import LT.box as B

import sys


#%% setup the PYTHON path

def add_sys_path(new_path):
    """ AddSysPath(new_path): adds a directory to Python's sys.path

    Does not add the directory if it does not exist or if it's already on
    sys.path. Returns 1 if OK, -1 if new_path does not exist, 0 if it was
    already on sys.path.
    """
    import sys, os

    # Avoid adding nonexistent paths
    if not os.path.exists(new_path): return -1

    # Standardize the path. Windows is case-insensitive, so lowercase
    # for definiteness.
    new_path = os.path.abspath(new_path)
    if sys.platform == 'win32':
        new_path = new_path.lower(  )

    # Check against all currently available paths
    for x in sys.path:
        x = os.path.abspath(x)
        if sys.platform == 'win32':
            x = x.lower(  )
        if new_path in (x, x + os.sep):
            return 0
    sys.path.append(new_path)
    return 1

orbit_mod90_python = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/python_modules'

if (add_sys_path(orbit_mod90_python)) < 0 :
    print(f'Cannot add {orbit_mod90_python} to PYTHONPATH as it does not exist !')
    sys.exit(-1)

#%% load the bundle stuff

import load_bundle as LB

b0 = LB.trajectory_bundle('Detector1.npz')

#%% get B-field in cart. corrds as the detector location

Bf = b0.get_B_fields(0, cartesian=True)

# B-field at detector in std coord. system
B_det = Bf[:,0]

#%% transform to local 
