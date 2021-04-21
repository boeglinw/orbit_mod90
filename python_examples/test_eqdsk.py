import numpy as np
import LT.box as B
import copy as C

import fluxpy as FL
# efit_directory = '../../orbit_public/MAST_efit/'
# efit_gfile_name = '029881.00252.dat'
#%%
efit_directory = '../../MAST-U/efit'
# efit_gfile_name = '29904_0.25_gfile-from-mast_mastu_script_div_nfil_manual-lim.eqdsk'
efit_gfile_name = '29904_0.25.eqdsk'
ret0 = FL.flux.get_flux(efit_directory, efit_gfile_name)

sel = FL.flux_par_mod.rlim != 0.

rlim0 = C.copy(FL.flux_par_mod.rlim[:FL.flux_par_mod.limitr])
zlim0 = C.copy(FL.flux_par_mod.zlim[:FL.flux_par_mod.limitr])

# plot the limiter
B.pl.plot(rlim0,zlim0)


#%%
efit_directory = '../../orbit_public/NSTX-U_efit/'
# efit_gfile_name = '29904_0.25_gfile-from-mast_mastu_script_div_nfil_manual-lim.eqdsk'
efit_gfile_name = '13511.004000_5.5kG_800kA_A1.6_kappa2.7'
ret1 = FL.flux.get_flux(efit_directory, efit_gfile_name)
sel = FL.flux_par_mod.rlim != 0.

rlim1 = C.copy(FL.flux_par_mod.rlim[:FL.flux_par_mod.limitr])
zlim1 = C.copy(FL.flux_par_mod.zlim[:FL.flux_par_mod.limitr])

# plot the limiter
B.pl.plot(rlim1, zlim1)


