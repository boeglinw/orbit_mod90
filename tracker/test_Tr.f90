! example of using boris tracker module with efit file
program test_T
  use Tracker
  use helper_functions_mod

  implicit none
  
  ! file name

  character*80 :: t_name
  character*4 :: c_E
  
  real(kind = 8), dimension(3) :: r0, ur0, rt0
  real(kind = 8), dimension(3) :: v0, vt0, uv0, uvt0

  ! for timing purposes
  real(kind = 8) :: start_time, end_time, elapsed_time
  
  integer(kind = 4) :: i,k, n_calc, ntra

  ! for limiter
  character(len = 132) :: directory = '../example_data/'
  character(len = 132) :: filename = 'MASTLIMIT00.dat'

  !efit_gfile_name = '029881.00252.dat'  ! efit result used to calculate the example track
  !efit_directory = '../example_data/'

  efit_gfile_name = '029880.00234.dat'  ! efit result used to calculate the example track
  efit_directory = '../example_data/'

  
  ! initial position (corresponds to the track in the example directory)
  !r0(1) = 0.32623276891949599
  !r0(2) = 1.7257147211829555
  !r0(3) = 2.8862357467520010E-002

  r0(1) = 0.31166603363303275 
  r0(2) = 1.6379398472054623 
  r0(3) = 3.0013999999999999E-002

  ur0 = r0/vect_mag(r0)
  
  ! initial velocity unit vector
  ! uv0(1) = 0.43438672503661674
  ! uv0(2) = -0.6748680841524517  
  ! uv0(3) = 0.59654106489357639 

  uv0(1) = 0.38713118426212662
  uv0(2) = -0.64032130488612071  
  uv0(3) = 0.66341395273293668

  print *, 'cartesian vector lengths = ', vect_mag(r0), vect_mag(uv0)

  ! convert to toroidal coordinates for integrator in toroida coords
  rt0(1) = vect_mag(r0)
  rt0(2) = pol_angle(r0(1:2)) ! pass only elments 1 and 3 of r0
  rt0(3) = r0(3)

  print *, 'toroidal position = ', rt0
  
  uvt0(1) = ur0(1)*uv0(1) + ur0(2)*uv0(2)
  uvt0(2) = -ur0(2)*uv0(1) + ur0(1)*uv0(2)
  uvt0(3) = uv0(3)
  print *, 'toroidal velocity = ', uvt0
  vt0 = uvt0
  
  print *, 'toroidal vector lengths = ', sqrt(rt0(1)**2 + rt0(3)**2), vect_mag(uvt0)
  
  ! set important parameters
  particle_charge = 1.*ec
  ! particle_mass_amu = 1.007347
  particle_mass_amu = 1.00
  particle_energy_mev = 3.
  track_length = 10.
  step_size = 0.01
  time_reversed = .true.
  
  ! initialize flux, get the efit results
  call init_tracker()

  ! set initical velocity
  v0 = uv0 * vmag
  
  print*, 'Particle energy = ', particle_energy,' J ', particle_energy/mev,' MeV'
  print*, 'Particle velocity = ', vmag, ' m/s'
  print*, 'particle q/m = ', q_over_m
  print*, 'particle omega = ', omega
  print*, 'Number of steps = ', Nsteps
  print*, 'step size = ', step_size
  print*, 'time step size = ', time_step_tor
  
  call load_flux()
  if (.not. flux_ok) then
     print*, 'cannot load eqdsk '
     stop
  endif

  ! initialize the limiter
  call set_limiter_file_name(filename)
  call set_limiter_directory(directory)
  call init_limiter

  ! select the tracker
  selected_tracker = bulirsch_stoer_t 
  !t_name = 'track_g_Tr_BS.dat'
  t_name = 'track_g_Tr_BS_cart.dat'
  !selected_tracker = boris_t 
  !t_name = 'track_g_Tr_BO.dat'

  
  ! setlect to check the limiter
  check_limiter = .false.
  print_hit = .true.
  
  call show_tracker()
  
  ! calc a trajectory
  ! calculate ntraj trajectories
  ntra = 1
  call cpu_time(start_time)
  do k = 1, ntra
     n_calc =  get_trajectory(r0, v0)
  enddo
  call cpu_time(end_time)
  elapsed_time = end_time - start_time
  print *, 'last get_trajectory status = ', n_calc
  print *, 'time for ', ntra, ' trajectories = ', elapsed_time 
  print *, 'time/trajectory = ', elapsed_time/ntra 

  ! write the trajectory data into a file
  open(unit = 16, file = t_name, status = 'unknown')
  
  do i = 1, n_calc
     write(16,*) (trajectory(i,k), k = 1, 3)
  enddo
  close(16)
  
end program test_T
