module boris_mod
  use constants_and_masses_mod
  use helper_functions_mod
  use limiter_mod
  use control_mod
  ! boris tracker module
  ! this contains the b_push routine that advances the position and velocities for the
  ! requested number of steps

contains 

  function b_push( rs, vs, dt, nsteps, bfield, efield, track) result(n_calc)
    ! calculates trajectory and returns the actual number of steps calculated, this should
    ! be smaller than nsteps if a limiter is hit and check_limiter is set to .true.
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
    integer(kind = 4), intent(in) ::  Nsteps
    integer(kind = 4) :: n_calc ! calculated steps
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
    track(1,1:3) = r0
    track(1,4:6) = v0    
    do i = 1, Nsteps-1
       ! perform 1 Boris step
       vminus = v0 + q_over_m*dt/2.*Efield(r0)   ! q_over_m  from constants_and_masses_mod
       tvect = q_over_m*dt/2.*Bfield(r0)
       vprime = vminus + cross_product(vminus,tvect) ! cross_product from helper_functions_mod
       svect = 2.*tvect/(1. + dot_product(tvect, tvect))
       vplus = vminus + cross_product(vprime, svect)
       v1 = vplus + q_over_m*dt/2.*Efield(r0)
       r1 = r0 + v1*dt
       ! save position and velocity data
       track(i+1,1:3) = r1
       track(i+1,4:6) = v1
       if (check_limiter) then   ! if selected check a limiter hit
          if (hit_lim(r1)) then
             n_calc = i+1        ! adjust the total number of steps calculated
             return              ! if a limiter is hit return
          endif
       endif
       r0 = r1
       v0 = v1
    enddo
    n_calc = Nsteps
  end function b_push
  
  
end module boris_mod

