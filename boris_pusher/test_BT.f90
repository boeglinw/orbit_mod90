! example of using boris tracker module with efit file
program test_BT
  use Tracker

  implicit none
  
  ! file name

  character*80 :: t_name
  character*4 :: c_E
  
  real(kind = 8), dimension(3) :: r0
  real(kind = 8), dimension(3) :: v0, uv0
  integer(kind = 4) :: i,k, n_calc

  ! for limiter
  character(len = 132) :: directory = '../example_data/'
  character(len = 132) :: filename = 'MASTLIMIT00.dat'

  ! initial position (corresponds to the track in the example directory)
  r0(1) = 0.32623276891949599
  r0(2) = 1.7257147211829555
  r0(3) = 2.8862357467520010E-002
  
  ! initial velocity unit vector
  uv0(1) = 0.43438672503661674
  uv0(2) = -0.6748680841524517  
  uv0(3) = 0.59654106489357639 

  ! set important parameters
  particle_charge = 1.
  particle_mass_amu = 1.007347
  particle_energy_mev = 3.
  track_length = 5.
  step_size = 0.01
  time_reversed = .true.
  efit_gfile_name = '029881.00252.dat'  ! efit result used to calculate the example track
  efit_directory = '../example_data/'

  ! initialize flux, get the efit results
  call init_BT()

  print*, 'Particle energy = ', particle_energy,' J ', particle_energy/mev,' MeV'
  print*, 'Particle velocity = ', vmag, ' m/s'
  
  call load_flux()
  if (.not. flux_ok) then
     print*, 'cannot load eqdsk '
     stop
  endif

  ! initialize the limiter
  call set_limiter_file_name(filename)
  call set_limiter_directory(directory)
  call init_limiter

  ! setlect to check the limiter
  check_limiter = .true.

  ! set the initial velocity
  v0 = vmag*uv0
  ! calc a trajectory
  n_calc =  get_trajectory(r0, v0)
  t_name = 'track_g_BT.dat'

  
  print *, 'get_trajectory status = ', n_calc

  ! write the trajectory data into a file
  open(unit = 16, file = t_name, status = 'unknown')
  
  do i = 1, n_calc
     write(16,*) (trajectory(i,k), k = 1, 3)
  enddo
  close(16)
  
end program test_BT
