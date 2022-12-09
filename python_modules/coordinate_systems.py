#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 09:45:29 2020

Coordinate Systems Setups

coordinate systems:
    
    Tokamak Global Coordinate System:  standard coord. system
    ---------------------------------------------------------
        Cartesian coordinates:
    
        z-axis : perp. to mid plane
        x-axis : in mid-plane (phi = 0)
        y-axis : in mid-plane (phi = pi/2.)
    
        Spherical coordinates:
            
        phi    : azimuthal angle measure from x-axis in ccw direction
        theta  : polar angle with respect to the z-axis
    
    
        toroidal coordinates
        
        phi : also the toroidale angle 
        R   : radial position : distance from z-axis
        Z   : same as z in cart. coord
        
        theta_p : poloidal angle: angle betweem x-axis and r when projected into a poloidal plane
        
        r : radial distance from magnetic axis in poloidal plane


    Local Coordinate System:
    ------------------------
        Basically the standard coordinate system rotated around the z-axis  
      
        Cartesian coordinates:
           
           z-axis : perp. to mid plane
           x-axis : in mid plane, pointing from center stack to detector detector or detector array position towards
           y-axis : in mid plane forming a RH coord system
           
        Polar coordinate system:
           
        used to describe detector orientation in local coordinate system:
               
          the detector plane is a plane perpendicular to the x-z plane containing the detector direction
           
          theta : angle of the detector plane with z-axis
          phi   : angle between detector direction and the intersection of the detector plane and the x-z plane.
                  For theta = 0 this corresponds to a rotation around the x-axis by -phi              
           
    Detector Coordinate System:
    ---------------------------    
        Cartesian coordinates:
            Origin location: particle detector center
            
            z-axis : pointing from detector center to collimator center
            x-axis : forms RHS coord. system with y- and z-axis 
            y-axis : forms RHS coord. system with y- and z-axis
            
        This is the coordinate system in which the initial values of trajectories are calculated by trajectory_bundles
                     

@author: boeglinw
"""

import numpy as np
import numpy.linalg as L


twopi = 2.*np.pi

dtr = np.pi/180.


#-----------------------------------------------------------------------------
# helper functions
#-----------------------------------------------------------------------------


def get_angle(x,y):
    # get phi for a vector
    angle = (np.arctan2(y,x) + twopi)%twopi
    return angle

def save_arccos(x,y):
    a = np.arccos( np.clip(x/y, -1.,1.) )
    return a


def get_unit_vector(r):
    # return unit vector in the direction of v
    return r/L.norm(r)

#-----------------------------------------------------------------------------
# coordinate conversions
#-----------------------------------------------------------------------------

def cart_to_spherical(r):
    x,y,z = r
    # convert cartesian to spherical
    rv = np.sqrt(x**2 + y**2 + z**2)
    phi = get_angle(x, y)
    theta = save_arccos(z, rv)
    return rv, phi, theta

def cart_to_rzphi(r):
    x,y,z = r
    # convert rzphi to spherical
    R = np.sqrt(x**2 + y**2)
    Phi = get_angle(x, y)
    Z = z
    return R, Z, Phi

#-----------------------------------------------------------------------------
# rotation matrices
#-----------------------------------------------------------------------------

def Rx(x):
    # rotation around x-axis by x:
    M = np.array([[1,0,0],[0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]])
    return M

def Ry(x):
    # rotation around y-axos by y
    M = np.array([[np.cos(x),0,np.sin(x)],[0, 1, 0], [-np.sin(x), 0 , np.cos(x)]])
    return M

def Rz(x):
    # rotation around z-axis by x
    M = np.array([[np.cos(x),-np.sin(x),0], [np.sin(x), np.cos(x), 0], [0, 0, 1]])
    return M



#-----------------------------------------------------------------------------
#   position and direction class
#-----------------------------------------------------------------------------
class pos_dir:
    def __init__(self, pos_type = 'x_y_z', position = None, dir_type = 'x_y_z', direction = None, alpha = 0., position_local = None):
        # list ot position types and keywprds
        self.coord_types = ['x_y_z', 
                            'r_phi_theta', 
                            'r_z_phi']
        # list of direction types and keywords
        self.dir_types = ['x_y_z', 
                         'phi_theta', 
                         'local_theta_phi']
        
        # check for valif position type
        if not (pos_type in self.coord_types):
            print('unknown position type ', pos_type, ' available values :', self.coord_types)
            return
        if not (dir_type in self.dir_types):
            print('unknown direction type ', dir_type, ' available values :', self.dir_types)
            return
        self.set_positions = {'x_y_z':self.set_pos_xyz,
                             'r_phi_theta': self.set_pos_r_phi_theta,
                             'r_z_phi':self.set_pos_r_z_phi }
        
        self.set_directions = {'x_y_z':self.set_dir_xyz,
                               'phi_theta': self.set_dir_theta_phi,
                               'local_theta_phi':self.set_dir_local_theta_phi}
        # default values
        self.pos_local = np.zeros(3)
        self.pos = np.array([1., 0., 0.])
        self.dir = np.array([-1., 0., 0.])
        self.R = 1.
        self.Phi = 0.
        
        self.alpha = alpha # detector head arm rotation
        # setup position if given
        if position is not None:
            self.set_positions[pos_type](*position)
        # setup direction if given
        if direction is not None:
            self.set_directions[dir_type](*direction)
        # setup position in local coord. system in given, this requried cartesian coordinates 
        if position_local is not None:
            self.set_pos_local(*position_local)
            
    #--------------------------------------------------------------------------
    # position functions
    #--------------------------------------------------------------------------
    def set_pos_r_z_phi(self, r, z, phi):
        self.Z = z
        self.R = r
        self.Phi = phi
        # calculate the position vector from r, z, phi
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        self.pos = np.array([x,y,z])
        
    def set_pos_xyz(self, x,y,z):
        self.pos = np.array([x,y,z])
        self.R = np.sqrt(x**2 + y**2)
        self.Z = z
        self.Phi = get_angle(x,y)
        

    def set_pos_r_phi_theta(self, r, phi, theta):
        print(r, phi, theta)
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        self.pos = np.array([x,y,z])
        self.R = np.sqrt(x**2 + y**2)
        self.Z = z
        self.Phi = get_angle(x,y)

    def set_pos_local(self, x,y,z):
        # set position of detector in local coordinate system
        self.pos_local = np.array([x,y,z])

    #--------------------------------------------------------------------------
    # direction functions (for detectors)
    #--------------------------------------------------------------------------
    def set_dir_theta_phi (self, theta, phi):
        # set direction using spherical coordinates
        self.theta = theta
        self.phi = phi 
        x = np.sin(theta)*np.cos(phi)
        y = np.sin(theta)*np.sin(phi)
        z = np.cos(theta)
        self.dir= np.array([x,y,z])
        
    def set_dir_xyz(self, x, y, z):
        # set direction along x,y,z calculating the unit vector along this direction
        self.dir = get_unit_vector(np.array([x,y,z]))

    #-------------------------------------------------------------------------
    # directions in local and detector system
    #-------------------------------------------------------------------------
    def set_dir_local_theta_phi(self, theta, phi):
        # set the direction relative to the local coordinate system
        # to use this one needs the orientation of the local coordinate system
        # x in the same direction as r-direction
        # y parallel to mid-plane and perp ccw to x
        # z like regular z-direction
        # phi rotation angle on plane that makes an angle theta wrsp to z- axis 
        # theta = 0  : plane is parallel to y-z plane
        # theta = 90 : plane is parallel to x-y plane
        #
        # alpha : final rotation around x-axis 
        self.theta_l = theta
        self.phi_l = phi
        x_l = -np.cos(phi)*np.sin(theta)
        y_l = np.sin(phi)
        z_l = np.cos(phi)*np.cos(theta)
        self.dir_local = np.array([x_l, y_l, z_l])
        # calculate
        
    def set_dir_local_xyz(self, x, y, z):
        # set direction in catesian corrdinates
        # create unit vector
        self.dir_local = get_unit_vector(np.array([x,y,z]))
        r_xz = np.sqrt(x**2 + z**2)
        self.theta_l = save_arccos(z, r_xz)
        # calculate phi
        # rotate into yz plane
        dvz = np.dot(Ry(self.theta_l), self.dir_local)
        self.phi_l = get_angle(dvz[1], dvz[2])
        
    #--------------------------------------------------------------------------
    # setup transformation matrices
    #--------------------------------------------------------------------------
    def set_detector_to_local(self):
        # transformation from detector coord system to local coord. system
        self.T_det_to_local_no_rotation = np.dot(Ry(-self.theta_l), Rx(-self.phi_l))
        self.T_arm_rotation = Rx(self.alpha)
        self.T_det_to_local = np.dot(self.T_arm_rotation, self.T_det_to_local_no_rotation)
        
    def set_detector_to_std(self):
        # transformation from detector coord. system to standard
        # make sure matrices are updated
        self.set_detector_to_local()
        self.set_local_to_standard()
        self.T_det_to_std = np.dot(self.T_local_to_std, self.T_det_to_local)
    
    def set_local_to_standard(self):
        # transformation from local (RP arm) coord. system to standard.
        # rotation around z-axis by toroidal angle Phi
        self.T_local_to_std = Rz(self.Phi)
        self.T_local_to_std_with_rotation = np.dot(self.T_local_to_std, self.T_arm_rotation)
        
    #--------------------------------------------------------------------------
    # vector transformation functions (rotations only)
    #--------------------------------------------------------------------------
    def detector_to_std(self, v):
        # given vector v in detector coord. syste, return vector in standard system
        return np.dot(self.T_det_to_std, v)
    
    def detector_to_local(self, v):
        # given vector v in detector coord. syste, return vector in local system
        return np.dot(self.T_det_to_local, v)

    def detector_to_local_no_rotation(self, v):
        # given vector v in detector coord. syste, return vector in local system, ingnoring arm rotation
        return np.dot(self.T_det_to_local_no_rotation, v)

    def local_to_std_with_rotation(self, v):
        # given vector v in local coord. syste, return vector std system including arm rotation
        return np.dot(self.T_local_to_std_with_rotation, v)
    
    def local_to_standard(self, v):
        # given vector v in detector coord. syste, return vector in standard system
        return np.dot(self.T_local_to_std, v)
        
        