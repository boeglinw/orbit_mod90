module em_fields_mod
  use orbit_parameters_mod
  use flux_par_mod
  use helper_functions_mod
  use control_mod
  implicit none

contains

  function bfield (r_loc, z_loc)  result(bf)
    !**********************************************************************
    !
    !     SUBPROGRAM DESCRIPTION:
    !          bfield returns the interpolated field
    !
    !     CALLING ARGUMENTS:
    !           r_loc.........r-position (input)
    !           z_loc.........z-position (input)
    !
    !     RETURNS
    !           bf( bpolr, bpolz, bphi, btotal, psirel) (output)
    !
    !           bpolr........r-Bfield component
    !           bpolz........z-Bfield component
    !           bphi.........phi-Bfield component
    !           btotal........Bfield magnitude
    !
    !          09/22/2018........changed top f90 W. Boeglin
    !          06/03/2020........rewritten W. Boeglin

    !
    !**********************************************************************


    implicit none

    ! position where field is required
    real(kind = 8), intent(in) :: r_loc, z_loc

    ! magnetic field in toroidal coord. system
    real(kind = 8), dimension(5)  :: bf


    ! interpolation function
    real(kind = 8):: seval
    ! real(kind = 8):: spline_eval

    ! local variables
    logical:: debug = .false.
    logical:: info = .false.

    logical:: is_outside_grid, is_inside, in_closed_flux_surface

    ! field components
    real(kind = 8) :: bpolr, bpolz, bphi
    ! magnitude
    real(kind = 8) :: btotal

    real(kind = 8) :: psi_rel
    real(kind = 8) :: bpol, fpnow

    ! local value of flux and its derivatives (from interpolation)
    real(kind = 8):: psi_loc, dpsi_dr, dpsi_dz


    integer(kind = 4):: ier, limfag, zeronow

    ! interpolation return array
    real(kind=8), dimension(6):: pds

    if (debug) print *, ' Get field at  R = ', r_loc, ' Z = ', z_loc

    ! check if point is inside grid from EQDSK file
    is_outside_grid = (r_loc.lt.rgrid(1)) .or. (r_loc.gt.rgrid(mw)) .or. (z_loc.lt.zgrid(1)) .or.  (z_loc.gt.zgrid(mh))

    if ( is_outside_grid ) then
       ! come here if point is out of field grid
       if (info) then
          print *, ' In flux, (R,Z) is outside grid: (',r_loc, z_loc,'), set field to 1e-8'
          print *, ' Grid R-limits : ', rgrid(1), rgrid(mw)
          print*,  ' Grid Z-limits : ', zgrid(1), zgrid(mh)
       endif
       bpolr = 1e-8
       bpolz = 1e-8
       bphi = 1.e-8
       btotal = 2.e-8

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
    dpsi_dr = pds(2)
    dpsi_dz = pds(3)

    ! poloidal field
    bpolz = dpsi_dr / r_loc
    bpolr = - dpsi_dz / r_loc
    bpol = sqrt (bpolr * bpolr + bpolz * bpolz)


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

    if (abs (current) .gt.10.) then
       if ( in_closed_flux_surface .and. (psi_rel .le. 1.) ) then
          !     take this branch if inside limiter and psi small enough i.e. inside plasma boundary
          fpnow = seval(mwfpol, psi_rel, xxxsi, fpol, bfpol, cfpol, dfpol)
          ! fpnow = spline_eval(mwfpol, psi_rel, xxxsi, fpol, bfpol, cfpol, dfpol)
       else
          fpnow = fpol (mwfpol)
       endif
    else
       fpnow = rzero * bzero
    endif

    ! toroidal field
    bphi = fpnow / r_loc

    !----------------------------------------------------------------------
    !   Total B field
    !----------------------------------------------------------------------

    btotal = sqrt (bphi * bphi + bpol * bpol)

    bf(1) = bpolr
    bf(2) = bpolz
    bf(3) = bphi
    bf(4) = btotal
    bf(5) = psi_rel

    if (time_reversed) then ! reverse field for time reversed calculations
       bf(1:3) = -bf(1:3)
    endif

    return


  end function bfield


  function bfield3(r) result(b)
    ! get the b-filed from the equilibrium in cartesian corrdinates at position r = (x,y,z)
    !
    ! this usese the function bfield
    implicit none

    real(kind=8), dimension(3), intent(in) :: r
    real(kind=8), dimension(3) :: b

    real(kind = 8), dimension(5) :: b_vect
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
    b_loc = b_vect(1)*u_r + b_vect(2)*u_z + b_vect(3)*u_phi

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

  function get_psi (r_loc, z_loc)  result(psi_rel)
      !**********************************************************************
      !
      !     SUBPROGRAM DESCRIPTION:
      !          get_psi returns relative flux at r_loc,z_loc
      !
      !     CALLING ARGUMENTS:
      !           r_loc.........r-position (input)
      !           z_loc.........z-position (input)
      !
      !     RETURNS
      !           psi_rel
      !
      !     psi_rel is defined as  psi_rel = ( psi(r,z) - psi(mag. axis) ) / ( psi(last closed flux surface) - psi(mag. axis) )
      !
      !     on axis : psi_rel = 0
      !     on last closed surface : psi_rel = 1.
      !
      !     0 < psi_rel < 1.
      !
      !**********************************************************************
  
  
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

  function get_psi_array(r, z, n) result(psi_a)

    ! calculate the relative flux for an array of points
    implicit none

    integer(kind = 4), intent(in) :: n
    ! position where field is required
    real(kind = 8), dimension(0:n-1), intent(in) :: r, z    

    ! results
    real(kind = 8), dimension(0:n-1) :: psi_a

    !  local variables
    integer(kind = 4) :: i

    do i = 0, n-1
       psi_a(i) = get_psi (r(i), z(i))
    enddo

    return

  end function get_psi_array
  
  
end module em_fields_mod
