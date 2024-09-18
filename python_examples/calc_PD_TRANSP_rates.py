#
# calculate the rates for  PD using trajectory bundles calculated with
# orbit_mod90 and calc_orbits.py (and derivatives)
#
import numpy as np
import os
import sys

import scipy.interpolate as SI

import orbit_view_data as vd

from LT import parameterfile as PF

import argparse as AG


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
#-------------------------------------------------------------

#-------------------------------------------------------------
# information class for TRANSP data
#-------------------------------------------------------------
class TRANSP_data:
    def __init__(self):
        self.dir = './TRANSP/'
        self.BT_file_name = None
        self.BB_file_name = None
        self.BB_var_type = 'g3_BB'
        self.BT_var_type = 'g3_BT'
    def show(self):
        print('dir = ', self.dir)
        print('BT_file_name = ', self.BT_file_name) 
        print('BB_file_name = ', self.BB_file_name) 
        print('BB_var_type = ', self.BB_var_type)
        print('BT_var_type = ', self.BT_var_type) 


#----------------------------------------------------------------------
# TRANSP use transp output
#----------------------------------------------------------------------
def EM_transp_init(info):
    # rr,zr, yntt, ynbt, ynbb  = EM_transp_init(info)
    # load transp Em as a function of R and Z
    # f is the file name
    # regular grid file extracted from CDF TRANSP output
    # see example: load_TRANSP_data.py
    cm2m = 1.e-2
    conv = 1./(4.*np.pi*cm2m**3)
    BT_data = np.load(info.dir + info.BT_file_name)
    BB_data = np.load(info.dir + info.BB_file_name)
    #store the grid data
    rg = BT_data['rg']
    zg = BT_data['zg']
    BT= BT_data[info.BT_var_type]*conv
    BB= BB_data[info.BB_var_type]*conv
    NY = BT+BB
    rt = rg*cm2m
    zt = zg*cm2m
    # assume emissivity is into 4pi, the accepance is m**2 Sr so need to divide Em by 4pi
    # setup 3d interpolation (default cubic interpolation)
    # r-range
    rr = rt[:,0]
    # z-range
    zr = zt[0,:]
    f_N  = SI.RectBivariateSpline(rr,zr,NY)
    f_BT = SI.RectBivariateSpline(rr,zr,BT)
    f_BB = SI.RectBivariateSpline(rr,zr,BB)
    return f_N, f_BB, f_BT
    


r_maxis = 1.0
z_maxis = 0.


#----------------------------------------------------------------------
# read parameter file

parser = AG.ArgumentParser()
parser.add_argument("control_file", nargs = '?', help="Control file ", default = 'control.data')
args = parser.parse_args()


# open control file
c_file = args.control_file
cd = PF.pfile(c_file)

# read control file:
    

# output setup
result_dir = cd['results_dir']
result_file = cd['results_file']

check_dir (result_dir)


# TRANSP result directory
TRANSP_dir = cd['TRANSP_dir']
TRANSP_BB = cd['TRANSP_BB_file']
TRANSP_BT = cd['TRANSP_BT_file']


# orbit directory
orbit_dir = cd['orbit_dir']


view_names = [f.strip() for f in cd['views'].split(',')]

#--------------------------------------------------------------------
# load the PD_views
PD_views = [ vd.view(orbit_dir + f) for f in view_names]


# load TRANSP calc
info = TRANSP_data()
info.dir = TRANSP_dir

info.BB_file_name = TRANSP_BB
info.BT_file_name = TRANSP_BT


# initialize transp data
f_N, f_BB, f_BT = EM_transp_init(info)

# set the source function
def EM_transp(x, r, z):
    # select total neutron yield
    return f_N.ev(r,z)

#%% calculate rate
rates = []
for PD_v in PD_views:
    rates.append(np.sum([R for R in PD_v.rates(EM_transp, in_plasma = True)]))
rates = np.array(rates)


#%% write results
o = open(result_dir + result_file, 'w')

o.write('# TRANSP rates\n')
o.write('#! name[i,0]/ rate[f,1]/ \n')
for i,pv in enumerate(PD_views):
    o.write(f'{pv.d.filename} {rates[i]}\n')
o.close()

