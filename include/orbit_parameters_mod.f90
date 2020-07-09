! module for orbit parameters and arrays
module orbit_parameters_mod
  
  implicit none
  ! parameters for R,Z grid
  integer(kind = 4), parameter :: nw=300
  integer(kind = 4), parameter :: nh=300
  
  integer(kind = 4), parameter :: nwnh=nw*nh
  integer(kind = 4), parameter :: nwrk=2*(nw+1)*nh
  
  ! parameters for boundary and limiter data 
  integer(kind = 4), parameter :: mbdry=1500
  integer(kind = 4), parameter :: nlimit=150
  
  ! parameters for spline interpolation
  integer(kind = 4), parameter :: kubicx = 4
  integer(kind = 4), parameter :: kubicy = 4
  integer(kind = 4), parameter :: lubicx = nw - kubicx + 1
  integer(kind = 4), parameter :: lubicy = nh - kubicy + 1
  
  integer(kind = 4), parameter :: nxknot = lubicx+2*kubicx-1
  integer(kind = 4), parameter :: nyknot = lubicy+2*kubicy-1
  
  ! max. number of detectors
  integer(kind = 4), parameter :: n_par = 20
  
  ! some control parameters
  integer(kind = 4), parameter :: n111 = 1
  integer(kind = 4), parameter :: n333 = 3
  
  ! io unit numbers
  integer(kind = 4), parameter:: neqdsk = 38

end module orbit_parameters_mod
