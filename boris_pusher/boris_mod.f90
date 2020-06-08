module boris_mod
  ! boris tracker module
  ! this contains the b_push routine that advances the position and velocities for the
  ! requested number of steps
  implicit none
  ! particle masses in kg
  real(kind = 8), parameter :: me = 9.10938356e-31
  real(kind = 8), parameter :: mp = 1.6726219e-27
  real(kind = 8), parameter :: mn = 1.674927471e-27
  real(kind = 8), parameter :: md = 3.343583772e-27
  real(kind = 8), parameter :: mamu = 1.66053892173e-27
  real(kind = 8), parameter :: ec = 1.60217662e-19
  real(kind = 8), parameter :: em_ratio = ec/me  ! em ratio in Kg/C
  real(kind = 8), parameter :: pi = 3.141592653589793
  real(kind = 8), parameter :: dtr = pi/180.
  real(kind = 8), parameter :: mev = ec*1.e6   ! MeV in J
  
  real(kind = 8) :: q  = 1.
  real(kind = 8) :: m = 1836.152
  real(kind = 8) :: q_over_m
  logical :: initialized = .false. 

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

  subroutine boris_init( qv,  mv )

    ! exmple : call init(1., 1836.152) 
    ! qv charge in e
    ! mv mass in electron masses e.g. proton
    
    real(kind = 8), intent(in) :: qv
    real(kind = 8), intent(in) :: mv
    q = qv
    m = mv
    q_over_m = q/m*em_ratio
    initialized = .true.
    print *, 'initialized boris with q = ', q, '(e) and m = ', m, ' (me),  ' , q_over_m, 'q/m ' 
    return
  end subroutine boris_init
  
  subroutine b_push( rs, vs, dt, nsteps, bfield, efield, track)  
    ! interfaces
    interface
       function Bfield(x) result(bf)
         real(kind = 8), dimension(3) :: bf
         real(kind = 8), dimension(3), intent(in) :: x
       end function Bfield
       
       function Efield(x) result(ef)
         real(kind = 8), dimension(3) :: ef
         real(kind = 8), dimension(3), intent(in) :: x
       end function Efield
    end  interface

    ! initial position
    real(kind = 8), dimension(3), intent(in):: rs
    ! intital velocity
    real(kind = 8), dimension(3), intent(in):: vs
    ! time step
    real(kind = 8) dt
    ! total number of steps
    integer(kind = 4) Nsteps
    ! track data
    real(kind = 8), dimension(nsteps,6), intent(out) :: track
    
    ! local variables
    ! loop variables
    integer(kind = 4) i, j, k

    real(kind = 8), dimension(3) :: r0
    real(kind = 8), dimension(3) :: r1
    real(kind = 8), dimension(3) :: v0
    real(kind = 8), dimension(3) :: v1
    real(kind = 8),dimension(3) :: vminus, vplus 
    real(kind = 8),dimension(3) :: tvect, svect
    real(kind = 8),dimension(3) :: vprime    

    ! save first position
    r0 = rs
    v0 = vs
    do i = 1, Nsteps
       ! perform 1 Boris step
       vminus = v0 + q_over_m*dt/2.*Efield(r0)
       tvect = q_over_m*dt/2.*Bfield(r0)
       vprime = vminus + cross_product(vminus,tvect)
       svect = 2.*tvect/(1. + dot_product(tvect, tvect))
       vplus = vminus + cross_product(vprime, svect)
       v1 = vplus + q_over_m*dt/2.*Efield(r0)
       r1 = r0 + v1*dt
       ! save position and velocity data
       do j = 1, 3
          track(i,j) = r0(j)
          track(i,j+3) = v0(j)
       enddo
       r0 = r1
       v0 = v1
    enddo
  end subroutine b_push
  
end module boris_mod

