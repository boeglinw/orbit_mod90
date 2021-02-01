program test_cylinder
  use boris_cylinder

  implicit none

  real(kind = 8), dimension(3) :: r_s
  real(kind = 8), dimension(3) :: v_s
  real(kind = 8), dimension(3) :: B_loc

  real(kind = 8) :: K, v, R_loc, d_loc, step_loc
  real(kind = 8) :: mev = 1.6e-13
  real(kind = 8) :: mp = 1.6727e-27
  integer(kind = 4) :: np
  
  K = 3.*mev
  v = sqrt(2.*K/mp)
  
  v_s = (/ 0.d0, 0.d0, v /)
  r_s = (/ -0.009d0,0.d0,0.d0 /)
  B_loc = (/ 0.d0, 0.3d0, .3d0 /)
  
  R_loc = 0.01d0
  d_loc = 0.05d0

  step_loc = 0.4d-3

  call set_b0(B_loc)
  
  call init(1.d0, 1836.152d0, R_loc, d_loc, step_loc, v_s )

  np = track_cylinder(r_s)  

  print *, 'calculated ', np, ' out of ', Nsteps, ' points'
 
end program test_cylinder
