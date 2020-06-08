! test boris with uniform field in z-direction
module Tracker
  use boris_mod
  implicit none
  real(kind = 8) :: particle_charge = 1.               ! particle charge in units of e
  real(kind = 8) :: particle_mass_amu  = mp/mamu       ! particle mass in amu default proton
  real(kind = 8) :: particle_energy_mev = 3.           ! particle energy in MeV
  real(kind = 8) :: particle_mass   = mp               ! particle mass in Kg default proton
  real(kind = 8) :: particle_energy =  3.*mev          ! particle energy in J
  real(kind = 8) :: particle_mass_me                   ! particle mass in me
  real(kind = 8) :: vmag                       ! particle velocity in m/s
 
  character(len = 132) :: efit_gfile_name = ''
  character(len = 132) :: efit_directory = ''
  logical :: time_reversed = .true.      ! if true invert B field for time reversed calculation
  logical :: flux_ok = .false.           ! if true B field calculation should be ready

  real(kind = 8) :: track_length = 10.   ! max. traj. length in m
  real(kind = 8) :: step_size = 0.01     ! step size in m
  real(kind = 8) :: time_step = 1e-3     ! time step should be step_size/vmat
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
    character(len = 132), intent(in) :: gfile_name
     efit_gfile_name = gfile_name
  end subroutine set_gfile_name

  subroutine set_directory(directory)
    implicit none
    character(len = 132), intent(in) :: directory
    efit_directory = directory
  end subroutine set_directory
  
  subroutine init_BT()

    implicit none
    ! set important values to the correct units
    particle_mass = particle_mass_amu*mamu        ! part. mass in Kg
    particle_energy = particle_energy_mev*mev     ! particle energy in J

    ! prepare input parameters for calculation
     ! calculate particle velocity
    vmag = sqrt(2.*particle_energy/particle_mass)

    ! total number of steps
    Nsteps = nint(track_length/step_size)

    ! particle mass in me
    particle_mass_me = particle_mass/me

    ! set values for tracker
    call boris_init(particle_charge, particle_mass_me)
    
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
    return
  end subroutine init_BT

  subroutine reset_BT()
    implicit none
    if (allocated(trajectory)) deallocate(trajectory)
    return
  end subroutine reset_BT
  
  function get_trajectory(r0, v0) result(ier)
    implicit none
    real(kind = 8), dimension(3) :: r0  ! initial position vector    
    real(kind = 8), dimension(3) :: v0  ! initial velocity vector
    integer(kind = 4) :: ier
        
    integer(kind = 4) :: i, j

    ier = 0
    if (.not. flux_ok) then
       print*, ' Flux not OK, cannot calculate B-field !'
       ier = -1
       return
    endif
    if (allocated (trajectory) ) then
       call b_push( r0, v0, time_step, Nsteps, bfield3, efield, trajectory)
       ier = 0
    else
       print *, 'trajectory not allocated cannot continue !'
       ier = -2
    endif
    return
  end function get_trajectory

  function bfield3(r) result(b)
    use boris_mod
    implicit none
    interface
       function bfield(rs,zs)
         real(kind = 8), dimension(4) :: bfield
         real(kind = 8), intent(in) :: rs
         real(kind = 8), intent(in) :: zs
       end function bfield
    end interface
    real(kind=8), dimension(3), intent(in) :: r
    real(kind=8), dimension(3) :: b
    
    real(kind = 8), dimension(4) :: b_vect
    real(kind = 8), dimension(3) :: u_r, u_phi, u_z, b_loc
    
    real(kind = 8) r_cyl, z_cyl
    
    ! unit vectoris of cyl.system
    u_z(1) = 0.
    u_z(2) = 0.
    u_z(3) = 1.
    ! unit vector in r-direction
    u_r(1:2) = r(1:2)/vect_mag(r(1:2))
    u_r(3) = 0.
    ! unit vectori in phi-direction
    u_phi(1) = -u_r(2)
    u_phi(2) = u_r(1)
    u_phi(3) = 0.
    
    ! calc toroidal coords from r,z
    r_cyl = sqrt(r(1)**2 + r(2)**2)
    z_cyl = r(3)
    
    b_vect = bfield(r_cyl,z_cyl)
    ! convert to cartesian coords
    if (time_reversed) then
       b_loc = -(b_vect(1)*u_r + b_vect(2)*u_z + b_vect(3)*u_phi)
    else
       b_loc = (b_vect(1)*u_r + b_vect(2)*u_z + b_vect(3)*u_phi)
    endif
    b(1) = b_loc(1)
    b(2) = b_loc(2)
    b(3) = b_loc(3)
    return
  end function bfield3
  
  
  function efield(r) result(e)
    implicit none
    real(kind=8), dimension(3), intent(in) :: r
    real(kind=8), dimension(3) :: e
    
    e(1) = 0.
    e(2) = 0.
    e(3) = 0.
    
    return
  end function efield
  
end module Tracker
