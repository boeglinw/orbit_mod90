! Main python interface module for BT

! include module sources to make sure they are available in python
include 'orbit_parameters_mod.f90'  ! make the orbit parameters available
include 'flux_par_mod.f90'          ! make the efit results available
include 'boris_mod.f90'             ! boris pusher module
include 'BT.f90'                    ! boris tracker modulem this is the one to
                                    ! work with
