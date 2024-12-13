! test boris with uniform field in z-direction
module Tracker
  use boris_mod  ! boris pusher module
  use bs_mod     ! Bulirsch-Stoer & RK4 integrator module 
  use em_fields_mod ! use the em_fields module to get the magnetic field
  use limiter_control_mod 
  use control_mod
  implicit none
  ! available trackers
  integer(kind = 4), parameter :: boris_t = 1
  integer(kind = 4), parameter :: bulirsch_stoer_t = 2
  integer(kind = 4) :: selected_tracker = 1 ! default selected tracker
  
  real(kind = 8) :: particle_charge_ec = 1.            ! particle charge in units of e
  real(kind = 8) :: particle_charge   = ec             ! particle charge in C
  real(kind = 8) :: particle_mass_amu  = mp/mamu       ! particle mass in amu default proton
  real(kind = 8) :: particle_energy_mev = 3.           ! particle energy in MeV
  real(kind = 8) :: particle_mass   = mp               ! particle mass in Kg default proton
  real(kind = 8) :: particle_energy =  3.*mev          ! particle energy in J
  real(kind = 8) :: particle_mass_me                   ! particle mass in me
  real(kind = 8) :: vmag                       ! particle velocity in m/s
 
  character(len = 132) :: efit_gfile_name = ''
  character(len = 132) :: efit_directory = ''
  
  logical :: flux_ok = .false.           ! if true B field calculation should be ready

  real(kind = 8) :: track_length = 10.   ! max. traj. length in m
  real(kind = 8) :: step_size = 0.01     ! step size in m
  real(kind = 8) :: time_step = 1e-3         ! time step should be step_size/vmat
  real(kind = 8) :: time_step_tor = 1e-2     ! time step should be step_size for scaled toroidal calculation
  integer(kind = 4) :: Nsteps ! number of steps
  
! trajectory
  real(kind = 8 ), allocatable, dimension(:,:) :: trajectory
  
contains

  subroutine load_flux()
    
    integer(kind = 4) :: get_flux
    integer(kind = 4) :: iret
    logical :: OK

    OK = .true.
    
    if (efit_directory .eq. '') then
       print*, '===> efit directory not set !'
       OK = .false.
    endif
    if (efit_gfile_name .eq. '') then
       print*, '===> efit file name not set !'
       OK = .false.
    endif
    
    if (.not. OK) then
       flux_ok = .false.
       return
    else
       iret = get_flux(efit_directory, efit_gfile_name)
    endif
    flux_ok = (iret .eq. 0)
    
    return

  end subroutine load_flux

  subroutine set_gfile_name(gfile_name)
    implicit none
    character(len = *), intent(in) :: gfile_name
     efit_gfile_name = gfile_name
  end subroutine set_gfile_name

  subroutine set_efit_directory(directory)
    implicit none
    character(len = *), intent(in) :: directory
    efit_directory = directory
  end subroutine set_efit_directory

  
  subroutine init_tracker()

    implicit none
    ! set important values to the correct units
    particle_mass = particle_mass_amu*mamu        ! part. mass in Kg
    particle_energy = particle_energy_mev*mev     ! particle energy in J

    ! prepare input parameters for calculation
     ! calculate particle velocity
    vmag = sqrt(2.*particle_energy/particle_mass)

    ! total number of steps
    Nsteps = nint(track_length/step_size)

    ! particle charge
    particle_charge = particle_charge_ec * ec    ! convert to C

    ! set q/m for tracker
    q_over_m = particle_charge/particle_mass

    ! set Omega used for scaling the EOM basically 
    ! 1/(r B) where r is the larmor radius and B the magntic field perp. to v
    ! all quantities are in SI units
    
    omega = particle_charge/sqrt(2.*particle_energy*particle_mass)

    ! allocate trajectory array
    if (.not. allocated(trajectory) ) then
       allocate( trajectory(NSteps, 6) ) 
    else
       ! reallocate trajectory as the number of steps could have changed
       deallocate(trajectory)
       allocate( trajectory(NSteps, 6) ) 
    endif
    ! set time step
    time_step = step_size/vmag
    time_step_tor = step_size
    return
  end subroutine init_tracker

  subroutine show_tracker
    print *, 'Tr---->tracker = ', selected_tracker
    print *, 'Tr---->time_reversed = ', time_reversed
    print *, 'Tr---->check_limiter = ', check_limiter
    print *, 'Tr---->print_hit = ', print_hit
  end subroutine show_tracker
  
  subroutine reset_tracker()
    implicit none
    if (allocated(trajectory)) deallocate(trajectory)
    return
  end subroutine reset_tracker
  
  function get_trajectory(r0, v0) result(n_calc)
    implicit none
    real(kind = 8), dimension(3) :: r0  ! initial position vector    
    real(kind = 8), dimension(3) :: v0  ! initial velocity vector
    integer(kind = 4) :: n_calc
        
    integer(kind = 4) :: i, j

    n_calc = 0
    if (.not. flux_ok) then
       print*, ' Flux not OK, cannot calculate B-field !'
       n_calc = -1
       return
    endif
    if (allocated (trajectory) ) then
       ! use the selected tracker
       trajectory = 0.
       select case(selected_tracker)
       case(boris_t)
          n_calc = b_push( r0, v0, time_step, nsteps, bfield3, efield, trajectory)
       case(bulirsch_stoer_t)
          n_calc = bs_push( r0, v0, time_step, nsteps, trajectory)
       case default
          print *, '------>no tracker savailable for selection = ', selected_tracker
          print *, '------>nothing to calculate !! '
          n_calc = -3
       end select
    else
       print *, 'trajectory not allocated cannot continue !'
       n_calc = -2
    endif
    return
  end function get_trajectory

  
end module Tracker
