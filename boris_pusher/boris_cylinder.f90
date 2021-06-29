module boris_cylinder
  use constants_and_masses_mod
  
  implicit none
  
  logical :: initialized = .false.
  logical :: hit = .false.


  real(kind = 8) :: m = mp/me
  real(kind = 8), dimension(3) :: B0 = 0.
  real(kind = 8), dimension(3) :: E0 = 0.  
  ! geometry
  
  real(kind = 8) :: C_r   ! Collimator radius
  real(kind = 8) :: R     ! cylinder radius
  real(kind = 8) :: d     ! length of hole
  real(kind = 8) :: step_size
  ! intital velocity
  real(kind = 8), dimension(3) :: vs
  
  logical :: stop_at_hits = .True. ! stop calculating the trajectory when hittin wall or collimator

  ! trajectories
  real(kind = 8 ), allocatable, dimension(:,:) :: track
  ! Number of steps
  integer(kind = 4) :: Nsteps
  
contains 

  function norm(v)
    real*8, dimension(:), intent(in) ::v
    real*8 norm
    norm = sqrt(dot_product(v,v))
  end function norm

  function cross_product(a, b) result(cp)
    real(kind = 8), dimension(3) :: cp
    real(kind = 8), dimension(3), intent(in) :: a, b
    
    cp(1) = a(2) * b(3) - a(3) * b(2)
    cp(2) = a(3) * b(1) - a(1) * b(3)
    cp(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  function get_par(a, b) result(a_par)
    ! calculate parallel component of vector a relative to vector b
    real(kind = 8), dimension(3), intent(in) :: a, b
    real(kind = 8), dimension(3) :: a_par
    a_par = b * dot_product(a,b)/norm(b)
  end function get_par

  function get_perp(a,b) result(a_per)
    ! calculate perpendicular component of vector a relative to vector b
    real(kind = 8), dimension(3), intent(in) :: a, b
    real(kind = 8), dimension(3) :: a_per
    a_per = a - get_par(a,b)
  end  function get_perp
  
  subroutine init( qv,  mv, Rv, Crv, dv, stepv, Nstepv, vsv )
    ! exmple : call init(1., 1836.152 ... 
    ! qv charge in e
    ! mv mass in electron masses e.g. proton
    ! Rv radius of cylinder
    ! Crv collimator radius
    ! dv distance collimator detector
    ! stepv size
    ! Nstepv number of steps
    ! vsv particle velocity
    
    real(kind = 8), intent(in) :: qv
    real(kind = 8), intent(in) :: mv
    real(kind = 8), intent(in) :: Rv
    real(kind = 8), intent(in) :: Crv
    real(kind = 8), intent(in) :: dv
    real(kind = 8), intent(in) :: stepv
    integer(kind = 4), intent(in) :: Nstepv
    ! intital velocity
    real(kind = 8), dimension(3), intent(in):: vsv

    real(kind = 8), dimension(3) :: v_par
    real(kind = 8), dimension(3) :: v_per

    ! total length
    real(kind = 8) :: l
    ! local variables


    real(kind = 8) :: w_c
    real(kind = 8) :: s
    real(kind = 8) :: R_c

    ! store input parameters
    q = qv  ! in ec
    m = mv  ! in me
    q_over_m = q/m*em_ratio

    R = Rv
    C_r = Crv
    
    d = dv
    step_size = stepv
    vs = vsv

    ! number of steps
    Nsteps = Nstepv

    ! only for non-zero B
    if (norm(B0) .gt. 1e-10) then
       ! cyclotron frequency
       w_c = q_over_M*norm(B0)
       ! velocity component parallel to B
       v_per = get_perp(vs, B0)
       ! cyclotron radius
       R_c = norm(v_per)/w_c

       ! estimate maximal path length
       s = sqrt(4.*R**2 + d**2)
       ! path length max.
       l = s*(1. + 1./24.*(s**2/R_c**2))
    else
       l = sqrt(4.*R**2 + d**2)
    endif
    

    ! allocate trajectory array
    if (.not. allocated(track) ) then
       allocate( track(NSteps, 6) ) 
    else
       ! reallocate trajectory as the number of steps could have changed
       deallocate(track)
       allocate( track(NSteps, 6) ) 
    endif

    initialized = .true.

    return
  end subroutine init
      
  function track_cylinder(rs) result(n_calc)
    ! special version of the boris tracker
    ! calculate the track starting at rs ( in the detector) until the collimator (a distance d away)
    ! check at every step if the track hits the wall
    ! initial position
    real(kind = 8), dimension(3), intent(in):: rs

    ! number of track data
    integer(kind = 4) :: n_calc ! calculated steps    
    ! loop variables
    integer(kind = 4) i, j, k

    real(kind = 8), dimension(3) :: r0
    real(kind = 8), dimension(3) :: r1
    real(kind = 8), dimension(3) :: v0
    real(kind = 8), dimension(3) :: v1
    real(kind = 8),dimension(3) :: vminus, vplus 
    real(kind = 8),dimension(3) :: tvect, svect
    real(kind = 8),dimension(3) :: vprime
 
    ! time step size
    real(kind = 8) :: dt
    ! distance from cylinder axis
    real(kind = 8) :: r_xy
    ! normalized time step
    real(kind = 8) :: qp    

    hit = .False.
    track = 0.

    ! save first position
    r0 = rs
    v0 = vs
    ! time step size
    dt = step_size/norm(v0)
    qp = q_over_m*dt/2.
    n_calc = 0

    ! if no field is present
    if ( (norm(E0) .lt. 1e-10) .and. (norm(B0) .lt. 1d-10 ) ) then
       ! no fields i.e. straight lines
       do i = 1, Nsteps
          r_xy = sqrt(r0(1)**2 + r0(2)**2)
          do j = 1, 3
             track(i,j) = r0(j)
             track(i,j+3) = v0(j)
          enddo
          if ( (.not. hit) .and. (r0(3) .le. d) .and. ((r_xy .gt. R) .or. (r_xy .gt. C_r)) ) then   ! cylinder wall or collimator was hit 
             hit = .True.
             if (stop_at_hits) then
                n_calc = i
                return
             endif
          endif
          r0 = rs + i*v0*dt 
       enddo
       n_calc = Nsteps
    else
    ! fields are present
       do i = 1, Nsteps
          ! save position and velocity data
          do j = 1, 3
             track(i,j) = r0(j)
             track(i,j+3) = v0(j)
          enddo
          r_xy = sqrt(r0(1)**2 + r0(2)**2)
          if ((.not. hit) .and. (r0(3) .le. d) .and. ((r_xy .gt. R) .or. (r_xy .gt. C_r))) then   ! cylinder wall or collimator was hit 
             hit = .True.
             if (stop_at_hits) then
                n_calc = i
                return
             endif
          endif
          ! perform 1 Boris step
          vminus = v0 + qp*E0
          tvect = qp*B0
          vprime = vminus + cross_product(vminus,tvect)
          svect = 2.*tvect/(1. + dot_product(tvect, tvect))
          vplus = vminus + cross_product(vprime, svect)
          v1 = vplus + qp*E0
          r1 = r0 + v1*dt
          r0 = r1
          v0 = v1
       enddo
       n_calc = Nsteps
    endif
  end function track_cylinder

  
  function B_const(r) result(B)
    real(kind = 8), dimension(3), intent(in) :: r 
    real(kind = 8), dimension(3) ::  B
    B = B0
  end function B_const
      
  function E_const(r) result(E)
    real(kind = 8), dimension(3), intent(in) :: r 
    real(kind = 8), dimension(3) :: E
    E = E0
  end function E_const

  subroutine set_B0(b)
    real(kind = 8), dimension(3), intent(in) :: b
    B0 = b
  end subroutine set_B0

  subroutine set_E0(e)
    real(kind = 8), dimension(3), intent(in) :: e
    E0 = e
  end subroutine set_E0
 
end module boris_cylinder

