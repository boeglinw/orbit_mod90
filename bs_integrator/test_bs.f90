program test_bs
  use bs_mod

  implicit none

  real(kind = 8), parameter :: me = 9.10938356e-31
  real(kind = 8), parameter :: mp = 1.6726219e-27
  real(kind = 8), parameter :: mn = 1.674927471e-27
  real(kind = 8), parameter :: md = 3.343583772e-27
  real(kind = 8), parameter :: mamu = 1.66053892173e-27
  real(kind = 8), parameter :: ec = 1.60217662e-19
  real(kind = 8), parameter :: em_ratio = ec/me  ! em ratio in Kg/C
  real(kind = 8), parameter :: mev = ec*1.e6   ! MeV in J

  real(kind = 8) :: eps = 1.d-5
  
  real(kind = 8) :: q = 1.*ec
  real(kind = 8) :: T = 3*mev
  real(kind = 8) :: v

  real(kind = 8), dimension(6) :: r, r0, rscale, dr
  real(kind = 8), dimension(6) :: drdt

  real(kind = 8) :: l = 1000.  ! path length
  real(kind = 8) :: dl = .005  ! step length
  real(kind = 8) :: l_tot, ds
  real(kind = 8) :: t0, dt, dt_done, dt_next ! time steps

  integer(kind = 4) :: nsteps, i, j, flag
  
  v = sqrt(2.*T/mp) ! proton velocity
  dt = dl/v

  nsteps = nint(l/dl)
   
  ! initial position
  r0(1:3) = 0.
  r0(4) = 0.
  r0(5) = 1.   ! initial direction in y-direction
  r0(6) = 0.
  t = 0.

  r = r0
  i = 0
  l_tot = 0.
  ! call lorentz(0.d0, r0, drdt, 6)
  ! print *, 'initial drdt = ', drdt
  do while ((l_tot .lt. l) .and. (i .lt. 100000))
     i = i + 1
     write(16,*) i, r(1:3), t, dt, l_tot 
     call bs_ode( t, dt, r, 6, lorentz, eps, flag)
     if (flag .ne. 0) then
        print *, ' flag = ', flag
     endif
     l_tot = l_tot + dt*v
  enddo
     
  stop 'all done'

  
  ! number of steps = 
  
contains

  function cross_product(a, b) result(cp)
    real(kind = 8), dimension(3) :: cp
    real(kind = 8), dimension(3), intent(in) :: a, b
    
    cp(1) = a(2) * b(3) - a(3) * b(2)
    cp(2) = a(3) * b(1) - a(1) * b(3)
    cp(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product
  
  subroutine lorentz(t, r, drdt, nn)
    ! generalized position r: x,y,z, vz, vy, vz
    ! use a b-field in the z-direction
    implicit none
    integer(kind = 4), intent(in) :: nn
    real(kind = 8), intent(in) :: t
    real(kind = 8), dimension(nn), intent(in) :: r
    real(kind = 8), dimension(nn), intent(out) :: drdt
    
    ! locals
    real(kind = 8), dimension(3) :: B = (/0.,0.,1./)  ! 1T in z-direction

    ! calculate derivatives from eq. of motion

    drdt(1) = r(4)*v
    drdt(2) = r(5)*v
    drdt(3) = r(6)*v
    drdt(4:6) = q/mp * cross_product(r(4:6), B)

    return
  end subroutine lorentz
    
end program test_bs
  
