module boris_tracker

  implicit none
  
  real(kind = 8) :: q  = 1.
  real(kind = 8) :: m = 1836.152
  real(kind = 8) :: q_over_m
  real(kind = 8), parameter :: em_ratio = 1.758820024e11  ! em ratio in Kg/C
  logical :: initialized = .false. 

contains 

  function cross_product(a, b) result(cp)
    real(kind = 8), dimension(3) :: cp
    real(kind = 8), dimension(3), intent(in) :: a, b
    
    cp(1) = a(2) * b(3) - a(3) * b(2)
    cp(2) = a(3) * b(1) - a(1) * b(3)
    cp(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product
  
  function norm(v)
    real*8, dimension(:), intent(in) ::v
    real*8 norm
    norm = sqrt(dot_product(v,v))
  end function norm

  subroutine init( qv,  mv )

    ! exmple : call init(1., 1836.152) 
    ! qv charge in e
    ! mv mass in electron masses e.g. proton
    
    real(kind = 8), intent(in) :: qv
    real(kind = 8), intent(in) :: mv
    q = qv
    m = mv
    q_over_m = q/m*em_ratio
    initialized = .true.

    return
  end subroutine init
  
  function b_push( rs, vs, l, dl, Bfield, Efield)  result(track)
    ! interfaces
    interface
       function cross_product(a,b) result(cp)
         real(kind = 8), dimension(3), intent(in) :: a
         real(kind = 8), dimension(3), intent(in) :: b
         end function cross_product
       function Bfield(x)
         real(kind = 8), dimension(3) :: Bfield
         real(kind = 8), dimension(3), intent(in) :: x
       end function Bfield
       
       function Efield(x)
         real(kind = 8), dimension(3) :: Efield
         real(kind = 8), dimension(3), intent(in) :: x
       end function Efield
    end  interface

    ! initial position
    real(kind = 8), dimension(3), intent(in):: rs
    ! intital velocity
    real(kind = 8), dimension(3), intent(in):: vs
    ! track data
    real(kind = 8), dimension(:,:), allocatable :: track
    
    ! total length
    real(kind = 8), intent(in):: l
    ! step size
    real(kind = 8), intent(in):: dl

    ! local variables
    integer(kind = 4) Nsteps
    ! loop variables
    integer(kind = 4) i, j, k

    real(kind = 8) dt
    real(kind = 8), dimension(3) :: r0
    real(kind = 8), dimension(3) :: r1
    real(kind = 8), dimension(3) :: v0
    real(kind = 8), dimension(3) :: v1
    real(kind = 8),dimension(3) :: vminus, vplus 
    real(kind = 8),dimension(3) :: tvect, svect
    real(kind = 8),dimension(3) :: vprime    

    ! now the code
    ! l and v need to have the same units
    dt = dl/norm(v0)
    ! number of steps
    Nsteps = nint(l/dl)

    ! allocate space for the track data
    allocate( track(Nsteps+1, 6) )

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
          track(i,j) = r0(i)
          track(i,j+3) = v0(i)
       enddo
       r0 = r1
       v0 = v1
    enddo
  end function b_push
  
end module boris_tracker

