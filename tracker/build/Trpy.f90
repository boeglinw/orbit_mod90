! Main python interface module for BT

! include module sources to make sure they are available in python
include 'orbit_parameters_mod.f90'       ! make the orbit parameters available
include 'flux_par_mod.f90'               ! make the efit results available
include 'limiter_control_mod.f90'        ! make the limiter controls available
include 'boris_mod.f90'                  ! boris pusher module
include 'bs_mod.f90'                     ! boris pusher module
include 'control_mod.f90'                ! control flags
include 'constants_and_masses_mod.f90'   !
include 'em_fields_mod.f90'              ! em fields routines
include 'Tr.f90'                         ! tracker modulem this is the one to
                                         ! work with
