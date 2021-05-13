#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:09:35 2020

acceptance calculations for reactangular and circular detector geometries

these are analytical/numerical (not based on Fibonacci grids)

@author: boeglinw
"""

import numpy as np
from scipy import integrate

#------------------------------------------------------------------------------
# helper functions
#------------------------------------------------------------------------------


def accept1d (xc, xd, d, xoff ):
    """
    Calculate 1d detecor acceptance

    Parameters
    ----------
    xc : Float
        collimator half width.
    xd : Float
        detector half width.
    d : Float
        distance coll. - detector.
    xoff : Float
        detector horizontal offset.

    Returns
    -------
    accept : Float
        integrated acceptance.

    """
    if (xc <= xd):
        full_opening = 2.*xc
        t_th1p = (xd - xc - xoff)/d
        t_th1n = (xc - xd - xoff)/d
    else:
        full_opening = 2.*xd
        t_th1n = (xd - xc - xoff)/d
        t_th1p = (xc - xd - xoff)/d
    t_th2p = (xc + xd - xoff)/d
    t_th2n = -(xc + xd + xoff)/d
    th1p = np.arctan(t_th1p)
    th1n = np.arctan(t_th1n)
    th2p = np.arctan(t_th2p)
    th2n = np.arctan(t_th2n)
    acc_x1 =  full_opening*(np.sin(th1p) - np.sin(th1n))
    acc_x2 =  (xc + xd + xoff)*(np.sin(th1n) - np.sin(th2n))
    acc_x3 = -d*(np.cos(th1n) - np.cos(th2n))
    acc_x4 =  (xc + xd - xoff)*(np.sin(th2p) - np.sin(th1p))
    acc_x5 =  d*(np.cos(th2p) - np.cos(th1p))
    accept =  acc_x1 + acc_x2 + acc_x3 + acc_x4 + acc_x5
    return accept

# special case for a hole i.e. the two radii as identical
def H_overlap(th, r,  d):
    """
    Calculate the overlap area for a hole of radius r and depth d

    Parameters
    ----------
    th : Float
        angle relative to center line.
    r : Float
        hole radius.
    d : Float
        length or depth of hole.

    Returns
    -------
    Float
        effective active area if this were a collimator that is seen by a detector for particles
        entering at an angle th relative to center line

    """
    x = d*np.tan(th)/(2.*r)
    if x >1 :
        return 0.
    T = 2./np.pi*(np.arccos(x) - x*np.sqrt(1-x**2))  # T is the transmission
    a = np.pi*r**2 * T
    return a

# more general case
# overlap of two circles R0 and R1 when shifted by d
def calc_overlap(R0, R1, d):
    """
    calculate the overlap of 2 circles of radius r1 and r2 separated by d:

    Parameters
    ----------
    R0 : Float
        radius of circle 1.
    R1 : Float
        radius of circle 2.
    d : Float
        distance between the centers.

    Returns
    -------
    Overlap area

    """
    r0 = max(R0,R1)  # r0 is the larger of the two
    r1 = min(R0,R1)  # r1 is the smaller of the two
    A1 = np.pi*r1**2
    if (d <= (r0 - r1)):
        # x0 = r1  # smaller is entirely inside the larger
        y0 = d
        return A1
    elif d > (r0 + r1):
        # x0 = 0.
        y0 = 0.
        return 0.
    # intersection point (right top)
    y0 = (r0**2 - r1**2 + d**2)/(2.*d)
    # x0 = np.sqrt(r0**2 - y0**2)
    #
    if d <= y0:
        x1_t = (y0 - d)/r1
        A1 = 0.5*r1**2*(2.*np.pi - 2.*np.arccos(x1_t) + 2.*x1_t*np.sqrt(1. - x1_t**2))
        x0_t = y0/r0
        A0 = 0.5*r0**2*(2.*np.arccos(x0_t) - 2.*x0_t*np.sqrt(1. - x0_t**2))
    elif d == y0:
        A1 = np.pi*r1**2/2.
        x0_t = y0/r0
        A0 = 0.5*r0**2*(2.*np.arccos(x0_t) - 2.*x0_t*np.sqrt(1. - x0_t**2))
    else:
        x1_t = (d - y0)/r1
        A1 = 0.5*r1**2*(2.*np.arccos(x1_t) - 2.*x1_t*np.sqrt(1. - x1_t**2))
        x0_t = y0/r0
        A0 = 0.5*r0**2*(2.*np.arccos(x0_t) - 2.*x0_t*np.sqrt(1. - x0_t**2))
    A = A0 + A1
    return A



#------------------------------------------------------------------------------
# detector geometry classes
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# general rectangular collimator - detector system
#------------------------------------------------------------------------------

class rect_detector:
    def __init__(self, coll_x = .1, coll_y = .1, det_x = .1, det_y = .1, dist = 1., off_x = 0., off_y = 0.):
        """
        Describes a rectangular collimator detector combination

        Parameters
        ----------
        coll_x : Float, optional
            half width of collimator in x-direction. The default is .1.
        coll_y : Float, optional
            half width of collimator in y-direction. The default is .1.
        det_x : Float, optional
            half width of detector in x-direction. The default is .1.
        det_y : Float, optional
            half width of detector in y-direction.. The default is .1.
        dist : Float, optional
            distance collimator-detector. The default is 1..
        off_x : Float, optional
            x-offset of detector wrsp to collimator. The default is 0..
        off_y : Float, optional
            y-offset of detector wrsp to collimator. The default is 0..

        Returns
        -------
        None.

        """
        # store detector parmeters
        self.c_x = coll_x # collimator half widths
        self.c_y = coll_y
        self.c_A = 4.*coll_x*coll_y
        self.d_x = det_x  # detector half width
        self.d_y = det_y
        self.d_A = 4.*det_x*det_y
        self.D = dist   # distance collimator - detector
        self.o_x = off_x # detector center offset relative to collimator center
        self.o_y = off_y
        # calculate the acceptance
        self.calc_acceptance()
        return

    def calc_acceptance(self):
        # calculate the acxep
        self.a_x = accept1d(self.c_x, self.d_x, self.D, self.o_x)
        self.a_y = accept1d(self.c_y, self.d_y, self.D, self.o_y)
        self.accept = self.a_x*self.a_y
# all done
#------------------------------------------------------------------------------
# cylindrical collimator detector- system
#------------------------------------------------------------------------------

class cyl_detector:
    def __init__(self, r = 1., d = 2.):
        """
        Describe a cylindrical collimator

        Parameters
        ----------
        r : TYPE, optional
            cylinder radius. The default is 1..
        d : TYPE, optional
            cylinder length. The default is 2..

        Returns
        -------
        None

        """
        self.r = r
        self.d = d
        self.th_max = np.arctan(2.*self.r/self.d)
        self.accept =  self.calc_acceptance()

    def calc_acceptance(self):
        self.f = lambda x: H_overlap(x, self.r, self.d)*np.sin(x)*np.cos(x)
        a_i = np.array(integrate.quad(self.f, 0., self.th_max)) * 2.*np.pi  # integrate of theta
        self.accept = a_i[0]
        self.daccept= a_i[1]


#------------------------------------------------------------------------------
# conical collimator - detector system
#------------------------------------------------------------------------------

class con_detector:
    def __init__(self, coll_r = .1, det_r = .1, dist = 1.):
        """
        calculate the acceptance of a conical collimator - detector system

        Parameters
        ----------
        coll_r : Float, optional
            collimator radius. The default is .1.
        det_r : Float, optional
            detector radius. The default is .1.
        dist : TYPE, optional
            distance collimator - detector. The default is 1..

        Returns
        -------
        None.

        """
        # store detector parmeters
        self.c_r = coll_r # collimator radius
        self.d_r = det_r  # detector radius
        self.D = dist   # distance collimator - detector
        self.th_max = np.arctan((self.c_r + self.d_r)/self.D)
        return

    def calc_acceptance(self):
        self.f = lambda x: calc_overlap(self.c_r, self.d_r, self.D*np.tan(x))*np.sin(x)*np.cos(x)
        a_i = np.array(integrate.quad(self.f, 0., self.th_max)) * 2.*np.pi
        self.accept = a_i[0]
        self.daccept = a_i[1]

