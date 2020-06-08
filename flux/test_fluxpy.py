import numpy as np
import LT.box as B

import fluxpy as FL
efit_directory = '../../orbit_public/MAST_efit/'
efit_gfile_name = '029881.00252.dat'

track_dir = '../../orbit_public/MAST_output/nml_orb_MAST_29881_252_g_ensemble/'
track_name = 'track_11111.data'

# get prev. caculated track
td = B.get_file(track_dir + track_name)
# get efit file
ret = FL.flux.get_flux(efit_directory, efit_gfile_name)

# tack locations
rt = td['r']
zt = td['z']

# fields
brt = td['br']
bzt = td['bz']
bphit = td['bphi']

# loop anc calculate points
#%%
B_calc = []
for i,rr in enumerate(rt):
    zz = zt[i]
    B_calc.append(FL.flux.bfield(rr,zz))
B_calc = np.array(B_calc)

br = B_calc[:,0]
bz = B_calc[:,1]
bphi = B_calc[:,2]
bmag = B_calc[:,3]


# test bfield_array
#%%
B_calc1 = FL.flux.bfield_array(rt, zt)
