#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 07:10:42 2021

fibonacci lattice class

e.g. http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
https://bduvenhage.me/geometry/2019/07/31/generating-equidistant-vectors.html

provide a lattice on a square, iin a circle or on a square

@author: boeglinw
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_points(ax, x,y,z, xlim = [-1.,1.], ylim = [-1.,1.], zlim = [-1., 1.]):

    ax.plot(x, y, z, '.')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)


class fibonacci_grid:
    def __init__(self, N = 20):
        #%% fibonaccig grid on unit square  
        self.GR = (1 + 5**0.5)/2.   # golden ratio
        self.N = N  # nu,mber of points
        self.calc_square()
    
    def set_npoints(self, N):
        self.N = N  # nu,mber of points
        self.calc_square()
        
    def calc_square(self):
        arg = np.arange(self.N)
        # points on unit square        
        self.uy = arg/self.N 
        self.ux = (arg/self.GR%1)
  
    def show_square(self):
        plt.figure()
        plt.plot(self.ux,self.uy,'.')
        plt.title(f'Unit square with {self.N} points')
        plt.xlabel('ux')
        plt.ylabel('uy')
        plt.gca().set_aspect('equal')
    
    def get_circle(self, R, N = None):
        if N is not None:
            self.set_npoints(N)
        # map grid onto a cicle
        self.theta_c = 2.*np.pi*self.ux
        self.r_c = np.sqrt(self.uy)
        return self.r_c*R, self.theta_c

    def show_circle(self):
        plt.figure()
        plt.plot( self.r_c*np.cos(self.theta_c), self.r_c*np.sin(self.theta_c), '.')
        plt.title(f'Unit circle with {self.N} points')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().set_aspect('equal')
        
    def get_sphere(self, theta_min, theta_max, N = None):
        # map onto a sphere
        # Lampert equal area transform
        if N is not None:
            self.set_npoints(N)
        self.phi_s = 2.*np.pi*self.ux
        delta_cth = np.cos(theta_min) - np.cos(theta_max)
        cth = np.cos(theta_min) - delta_cth*self.uy
        self.theta_s = np.arccos(cth)
        return self.phi_s, self.theta_s
    
    def show_sphere(self):    
        xss = np.cos(self.phi_s)*np.sin(self.theta_s)
        yss = np.sin(self.phi_s)*np.sin(self.theta_s)
        zss = np.cos(self.theta_s)
        ax1 = plt.figure().add_subplot(111, projection='3d')
        plot_points(ax1, xss, yss, zss)
        
# end of class