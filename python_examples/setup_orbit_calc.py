#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 13:15:18 2024

setup example orbit calculation in local directory


 >cp ../orbit_mod90/python_examples/calc_orbits.py ./python/.
 >cp ../orbit_mod90/python_examples/calc_orbit_control_4ch_29880.data ./data/.
 >cp ../orbit_mod90/python_examples/detector_head_4ch.data ./data/. 
 >cp ../orbit_mod90/example_data/g029880.00234.dat ./data/.
 >cp ../orbit_mod90/example_data/MASTLIMIT2013.dat ./data/.


@author: boeglinw
"""


import os
import sys
import shutil as SH


#-------------------------------------------------------------
# make directory if needed
def check_dir(out_dir):
    # create output directory (if necessary)
    if os.path.isdir(out_dir):
        print(70*'=')
        print(f'-----> {out_dir}  exists, will use it ')
        print(70*'=')
    else:
        try:
            print(70*'=')
            print(f'----> Try to create : {out_dir}')
            os.makedirs(out_dir, exist_ok=True) 
            print(70*'=')
        except Exception as msg:
            print("problem : ", msg)
            sys.exit()

#%% location of the orbit_mod90
orbit_mod90_dir = '/Users/boeglinw/Documents/boeglin.1/Fusion/Fusion_Products/orbit_mod90/'

example_dir = orbit_mod90_dir + '/python_examples/'
example_data = orbit_mod90_dir + '/example_data/'

#%% setup directories for running the examples

check_dir('./python/')
check_dir('./data/')


SH.copy(example_dir + 'calc_orbits.py', './python/.')

SH.copy(example_dir + 'calc_orbit_control_4ch_29880.data', './data/.')
SH.copy(example_dir + 'detector_head_4ch.data', './data/.')

SH.copy(example_data + 'g029880.00234.dat', './data/.')
SH.copy(example_data + 'MASTLIMIT2013.dat', './data/.')

#%%
print('all setup . Try to run the command : %run python/calc_orbits.py -c ./data/calc_orbit_control_4ch_29880.data')


