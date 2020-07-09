#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:41:45 2020

@author: boeglinw
"""

import numpy as np

f18 = np.loadtxt('fort.18')
f18o = np.loadtxt('../../orbit_public/temp/fort.18')


print('new data shape = ', f18.shape)

def compare(i, n):
    r = f18[:n,i]/f18o[:n,i]
    print ("max for i = " , i, ' = ', r.max() )
    
# all done

def show(i,n):
    for j, vv in enumerate(f18[:n]):
        print('i = ', i,' row = ',  j, ' values = ',  vv[i], f18o[j,i], ' r = ', vv[i]/f18o[j,i])
# all done
def show_one(i,j):
    print('i = ', i,' row = ',  j, ' values = ',  f18[j,i], f18o[j,i], ' r = ', f18[j,i]/f18o[j,i])
# all done
