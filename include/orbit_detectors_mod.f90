module orbit_detectors_mod
  use orbit_parameters_mod

  real(kind = 8), dimension(n_par):: theta_port
  real(kind = 8), dimension(n_par):: phi_port
  real(kind = 8), dimension(n_par):: gyro_angle
  real(kind = 8), dimension(n_par):: pitch_angle
  real(kind = 8), dimension(n_par):: rdist
  real(kind = 8), dimension(n_par):: ZDist
  real(kind = 8), dimension(n_par):: phdangle

  real(kind=8), dimension(n_par) :: ports
  real(kind=8), dimension(n_par) :: portsph
  real(kind=8), dimension(n_par) :: al0s
  real(kind=8), dimension(n_par) :: be0s
  real(kind=8), dimension(n_par) :: rds
  real(kind=8), dimension(n_par) :: zds
  real(kind=8), dimension(n_par) :: phda

  integer(kind = 4):: detec
  integer(kind = 4), dimension(n_par):: detector_number
  integer(kind = 4), dimension(n_par):: channel_number
  integer(kind = 4), dimension(n_par):: detector_id

end module orbit_detectors_mod
