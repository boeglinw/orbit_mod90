# integrate psi along the path
# use orbit output
# make sure that one obtains the same number if using the same function
# and the same parameters
#
# modification WB March 2015 include acceptance in rate calculation

import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as pl
# module to load trakectory bundles
import load_bundle as LB
#
# use the datafile version that can read the parameters
from LT.pdatafile import pdfile


class view:
    def __init__(self, file):
        # read data for the trajectory
        self.d = LB.trajectory_bundle(file)
    #end of init
    
    def rates(self, Em, in_plasma = True):
        # iterator to get the  rates for all trajectories
        counter = 0
        while counter < self.d.n_trajectories:
            t_data = self.d.get_trajectory(counter)
            x,y,z,r =  t_data[:4,:]
            br, bz, bpj, bt, psi_rel = self.d.get_B_fields(counter)
            # if selected use only the part of the trjectory inside the plasma
            if in_plasma:
                sl = (0. <= psi_rel) & (psi_rel <= 1.)
            else:
                sl = slice(0, x.size)
            pos_psi_rel = psi_rel[sl] 
            Sdl = Em(pos_psi_rel, r[sl], z[sl])*self.d.information['step_size']*self.d.acceptance[counter]
            Sdl_int = Sdl.sum()
            yield Sdl_int
            counter += 1        
