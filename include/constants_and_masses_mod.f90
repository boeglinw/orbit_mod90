module constants_and_masses_mod
  implicit none
  ! constants
  real(kind = 8), parameter :: pi  = 3.1415926535897932d0
  real(kind = 8), parameter :: twopi=2.0*pi
  real(kind = 8), parameter :: dtr = pi/180.
  
  ! particle masses in kg
  real(kind = 8), parameter :: me = 9.10938356e-31
  real(kind = 8), parameter :: mp = 1.6726219e-27
  real(kind = 8), parameter :: mn = 1.674927471e-27
  real(kind = 8), parameter :: md = 3.343583772e-27
  real(kind = 8), parameter :: mamu = 1.66053892173e-27
  ! elementary charge
  real(kind = 8), parameter :: ec = 1.60217662e-19
  real(kind = 8), parameter :: em_ratio = ec/me  ! em ratio in Kg/C
  ! energy
  real(kind = 8), parameter :: ev = ec   ! eV in J
  real(kind = 8), parameter :: mev = ec*1.e6   ! MeV in J

  ! conversion factor for magnetic moment
  real(kind = 8), parameter :: tmu = 2.0e-07
  
  ! default particle parameters
  real(kind = 8) :: q  = 1.*ec
  real(kind = 8) :: m_me = mp/me
  
  ! this varies depending on kinematics and needs to be adjusted
  real(kind = 8) :: q_over_m = ec/mp  ! default value for proton
  real(kind = 8) :: omega = 1.     ! omega is the particle gyrofrequency divided by the local magnetic
                                      ! field and the magnitude of the velocity
  
end module constants_and_masses_mod
