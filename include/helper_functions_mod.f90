module helper_functions_mod
  use constants_and_masses_mod
  implicit none
contains
  
  function cross_product(a, b) result(cp)
    real(kind = 8), dimension(3) :: cp
    real(kind = 8), dimension(3), intent(in) :: a, b
    
    cp(1) = a(2) * b(3) - a(3) * b(2)
    cp(2) = a(3) * b(1) - a(1) * b(3)
    cp(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product
  
  
  function vect_mag(v) result(v_mag)
    real*8, dimension(:), intent(in) ::v
    real*8 :: v_mag
    v_mag = sqrt(dot_product(v,v))
  end function vect_mag
  
  
  function pol_angle(rv) result(angle)
    ! calculate the polar angle for a 2d vector rv. The angle is between 0 and 2 pi
    real(kind = 8), dimension(2), intent(in) :: rv
    real(kind = 8) :: angle
    
    angle = mod(atan2(rv(2), rv(1) ) + twopi, twopi)
    return
  end function pol_angle
    
end module helper_functions_mod
  
