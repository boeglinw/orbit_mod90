
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
  !           bf( bpolr, bpolz, bphi, btotal) (output)
  !
  !           bpolr........r-Bfield component 
  !           bpolz........z-Bfield component 
  !           bphi.........phi-Bfield component
  !           btotal........Bfield magnitude 
  !                                                                  
  !     REFERENCES:                                                  
  !          (1)                                                     
  !          (2)                                                     
  !                                                                  
  !     RECORD OF MODIFICATION:                                      
  !          12/04/84..........first created                         
  !          21/11/86..........revised                               
  !          04/08/96..........revised by Q.Peng                     
  !          29/01/99..........revised by D. Darrow for FIGOL        
  !          09/22/2018........changed top f90 W. Boeglin
  !          06/03/2020........rewritten W. Boeglin

  !                                                                  
  !**********************************************************************
    
  use orbit_parameters_mod
  use flux_par_mod

  implicit none

  ! position where field is required
  real(kind = 8), intent(in) :: r_loc, z_loc

  ! magnetic field in toroidal coord. system
  real(kind = 8), dimension(4)  :: bf

  
  ! interpolation function
  real(kind = 8):: seval
  
  ! local variables
  logical:: debug = .false.
  logical:: info = .false.
  
  logical:: is_outside_grid, is_inside, in_closed_flux_surface

  ! field components
  real(kind = 8) :: bpolr, bpolz, bphi  
  ! magnitude
  real(kind = 8) :: btotal
  
  real(kind = 8) :: xsinow
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
  
  if (abs (sifb - sifm) .gt.1.e-8) then  
     xsinow = (psi_loc - sifm) / (sifb - sifm)  
  else  
     xsinow = 1.2       
  endif

  !     Check to see whether sample point is inside limiter region
  in_closed_flux_surface = is_inside (r_loc, z_loc, nbdry, rbdry, zbdry)
  ! check of plasma current (piii from eqdsk file) larger than 10A
  ! interpolate poloidal stream function

  if (debug) then
     print *, 'in close flux surface = ', in_closed_flux_surface
     print *, 'plasma current = ', current
     print *, 'relative flux = ', xsinow
  endif
  
  if (abs (current) .gt.10.) then       
     if ( in_closed_flux_surface .and. (xsinow .le. 1.) ) then  
        !     take this branch if inside limiter and psi small enough i.e. inside plasma boundary
        fpnow = seval(mwfpol, xsinow, xxxsi, fpol, bfpol, cfpol, dfpol)
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
  
  return  
  
  
end function bfield
