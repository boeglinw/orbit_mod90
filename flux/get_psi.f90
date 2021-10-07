
function get_psi (r_loc, z_loc)  result(psi_rel)
    !**********************************************************************
    !
    !     SUBPROGRAM DESCRIPTION:
    !          get_psi returns relative flux at point r_loc, z_loc
    !
    !     CALLING ARGUMENTS:
    !           r_loc.........r-position (input)
    !           z_loc.........z-position (input)
    !
    !     RETURNS
    !           psi_rel
    !**********************************************************************

    use orbit_parameters_mod
    use flux_par_mod
    use helper_functions_mod
    use control_mod

    implicit none

    ! position where field is required
    real(kind = 8), intent(in) :: r_loc, z_loc

    ! interpolation function
    real(kind = 8):: seval
    ! real(kind = 8):: spline_eval

    ! local variables
    logical:: debug = .false.
    logical:: info = .false.

    logical:: is_outside_grid, is_inside, in_closed_flux_surface


    real(kind = 8) :: psi_rel
    real(kind = 8) :: bpol, fpnow

    ! local value of flux and its derivatives (from interpolation)
    real(kind = 8):: psi_loc, dpsi_dr, dpsi_dz


    integer(kind = 4):: ier, limfag, zeronow

    ! interpolation return array
    real(kind=8), dimension(6):: pds

    if (debug) print *, ' Get psi_rel at  R = ', r_loc, ' Z = ', z_loc

    ! check if point is inside grid from EQDSK file
    is_outside_grid = (r_loc.lt.rgrid(1)) .or. (r_loc.gt.rgrid(mw)) .or. (z_loc.lt.zgrid(1)) .or.  (z_loc.gt.zgrid(mh))

    if ( is_outside_grid ) then
       ! come here if point is out of field grid
       if (info) then
          print *, ' In flux, (R,Z) is outside grid: (',r_loc, z_loc,'), set field to 1e-8'
          print *, ' Grid R-limits : ', rgrid(1), rgrid(mw)
          print*,  ' Grid Z-limits : ', zgrid(1), zgrid(mh)
       endif
       psi_rel  = 1e20  ! non-sens value for outside grid
       return
    else
       ! using n333 means icalc = 3, seva2d contains f, dfdx, dfdy returned in pds
       ! f    = pds(1)
       ! dfdx = pds(2)
       ! dfdy = pds(3)
       !
       ! bkx, lkx, bky, lky, c,  are interpolation parameters
       ! mw, mh the total grid size
       call seva2d (bkx, lkx, bky, lky, c, mw, mh, r_loc, z_loc, pds, ier, n333)
       ! call spline_eval_2d (bkx, lkx, bky, lky, c, mw, mh, r_loc, z_loc, pds, ier, n333)
    endif

    ! flux
    psi_loc = pds(1)

    ! calculate relative flux, interpolated value is stored in pds(1)
    ! sifb and sifm are from read_eqdsk

    if (abs (psi_bdry - psi_axis) .gt.1.e-8) then
       psi_rel = (psi_loc - psi_axis) / (psi_bdry - psi_axis)
    else
       psi_rel = 1.2
    endif

    !     Check to see whether sample point is inside limiter region
    in_closed_flux_surface = is_inside (r_loc, z_loc, nbdry, rbdry, zbdry)
    ! check of plasma current (piii from eqdsk file) larger than 10A
    ! interpolate poloidal stream function

    if (debug) then
       print *, 'in closed flux surface = ', in_closed_flux_surface
       print *, 'plasma current = ', current
       print *, 'relative flux = ', psi_rel
    endif

end function get_psi


