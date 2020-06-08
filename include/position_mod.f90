module position_mod
  ! contains the current position and velocity along a calculated orbit
  implicit none

  real(kind = 8), dimension(4):: r
  real(kind = 8), dimension(4):: v
  real(kind = 8), dimension(4):: b

  real(kind = 8) :: s
  real(kind = 8) :: sstp
  real(kind = 8) :: tol

end module position_mod

