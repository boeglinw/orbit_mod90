! Main python interface module

! include module sources to make sure they are available in python
include 'orbit_parameters_mod.f90'
include 'flux_par_mod.f90'

module flux
contains
  include 'get_flux.f90'
  include 'bfield.f90'
  include 'is_inside.f90'
  include 'checklim.f90'
end module flux

