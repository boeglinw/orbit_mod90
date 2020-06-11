program test_limiter

  use limiter_mod
  implicit none
  logical :: check
  real(kind = 8) :: r, z, phi
  
  limiter_directory = './'
  limiter_filename = 'MASTLIMIT00.dat'

  call init_limiter
  
  print *, ' enter r,z coords :'
  read (5,*) r, z, phi
  phi = phi*dtr
  
  check = limiter_hit(r, z, phi)
  print *, 'limiter_hit = ', check
  
end program test_limiter
