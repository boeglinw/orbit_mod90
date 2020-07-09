module limiter_mod
  use orbit_parameters_mod
  use constants_and_masses_mod
  use helper_functions_mod
  use control_mod

  implicit none

  
  !    ntorreg number of toroidal regions (ie intervals in phi, the toroidal
  !       angle).  must be .le. ntorrmx
  !     phi_start toroidal angle at which each toroidal region starts (ending
  !       angle is the same as the starting angle of the next region, and
  !       there is one extra element at the end of this array, which is 
  !       set, by definition, to 2*pi) !     ntorrmx is the maximum number of toroidal regions (ie intervals in
  !       phi, the toroidal angle).  used to dimension arrays.

  
  integer(kind = 4), parameter :: npolyvmx = 100 ! maximum number of points in a polygon
  integer(kind = 4), parameter :: npolymx = 10 ! max number pof polygons per toroidal region
  integer(kind = 4),parameter :: ntorrmx = 30  ! max number of toroidal regiuons
  
  type polygon
     integer(kind = 4) :: nv
     logical :: intersects_midplane
     real(kind = 8), dimension(npolyvmx):: x
     real(kind = 8), dimension(npolyvmx):: y
     real(kind = 8) :: x0_min !  min values for intersections where y = 0.
     real(kind = 8) :: x0_max !  max values for intersections where y = 0.
  end type polygon

  type tor_region
     real(kind = 8) :: phi_start         ! toroidal angle of region start
     integer(kind = 4) :: np             ! number of polygons
     type(polygon), dimension(npolymx) :: poly   ! array of polygons
  end type tor_region

  !     ntorreg: number of toroidal regions (ie intervals in phi, the toroidal
  !       angle).  must be .le. ntorrmx
  !     phi_start: toroidal angle at which each toroidal region starts (ending
  !       angle is the same as the starting angle of the next region, and
  !       there is one extra element at the end of this array, which is 
  !       set, by definition, to 2*pi)

  integer(kind = 4) :: ntorreg
  real(kind = 8), dimension(ntorrmx + 1) :: phi_start ! phy starting angles for toroidal regions
  type(tor_region), dimension(ntorrmx + 1) ::t_regions ! toroidal regions
  
  character(len = 132) :: limiter_filename = ''
  character(len = 132) :: limiter_directory = ''
  
  ! nlimf is the unit number for reading the limiter data file
  integer(kind = 4), parameter :: nlimf = 17

  logical :: limiter_initialized = .false.
  
contains


  
  function locate (z, zlist, nz) result(pos)

    !  This function searches through a list of z values (assumed
    !  sorted lowest to highest) to find the pair of z values in the list
    !  between which the given z lies.  The routine returns the index of
    !  the first of the pair of z's.  If the z value is out of range of
    !  the list, then an index of -1 is returned.

    implicit none
    
    real(kind = 8), intent (in) :: z
    real(kind = 8), dimension(nz), intent(in) :: zlist
    integer(kind = 4), intent(in) :: nz
    integer(kind = 4) pos

    ! locals
    integer (kind = 4) ::  top, bot, mid

    if ((z .gt. zlist(nz)) .or. (z .lt. zlist(1))) then
       pos = -1
       return
    endif
    
    !   set up binary search
    top = nz
    bot = 1
    mid = (nz + 1) / 2
    !    WHILE loop follows
    do while ((mid .ne. top) .and. (mid .ne. bot))
       if ((z .ge. zlist(mid)) .and. (z .le. zlist(mid+1))) then
          ! found the right interval, so return
          pos = mid
          return
       endif
       if (z .lt. zlist(mid)) then
          top = mid
       else             
          bot = mid
       endif
       mid = (top + bot) / 2
    enddo
    pos = mid
    if (mid .ge. nz) pos = nz - 1
    return
  end function locate

  function in_polygon(x, y, p) result(inside)
    ! function checks if point x, y is inside polygon p
    implicit none
    interface
       function is_inside (x, y, nd, x_p, y_p)  result(OK)
         integer, intent(in) :: nd
         real(kind = 8), dimension(nd), intent(in) :: x_p, y_p
         real(kind = 8), intent(in) :: x, y
         logical :: OK
       end function is_inside
    end interface

    real(kind = 8), intent(in) :: x
    real(kind = 8), intent(in) :: y
    type(polygon), intent(in) :: p
    logical :: inside

    inside = is_inside(x, y, p%nv, p%x, p%y)
    return
  end function in_polygon

  subroutine print_poly(c, p)
    implicit none
    type(polygon), intent(in):: p
    character(len = *), intent(in):: c

    integer(kind = 4) :: i

    print *,'---------------------------------------------------------------'
    print *,'comment = ', c
    print *, 'polygon number of vertices = ', p%nv
    do i = 1, p%nv
       print*, 'polygon point ', i, ' x = ', p%x(i), ',  y = ', p%y(i)
    enddo
    print *,'---------------------------------------------------------------'

  end subroutine print_poly
    
    
  subroutine get_mid_plane(p)
    ! find mid-plane intersection of polygon, is none exist the values for x0 are large negative numbers
    implicit none
    type(polygon), intent(inout):: p
    integer(kind = 4) :: i, j, nm
    integer(kind = 4) :: sign_old, sign_new
    real(kind = 8), dimension(:), allocatable :: rmid
    real(kind = 8) :: rmin, rmax
    real(kind = 8), dimension(2) :: x_loc, y_loc
    ! allocate space for temp .array
    allocate(rmid(p%nv))

    x_loc(1) = p%x(1)
    y_loc(1) = p%y(1)
    sign_old = int(sign(1.d0, y_loc(1)))
    nm = 0
    p%intersects_midplane = .false.
    do i = 2, p%nv + 1
       if (i .gt. p%nv) then   ! close polygon to get both limits
          x_loc(2) = p%x(1)    ! make the first vertex the new vertex after the last one
          y_loc(2) = p%y(1)
       else
          x_loc(2) = p%x(i)
          y_loc(2) = p%y(i)
       endif
       sign_new = int(sign(1.d0, y_loc(2)) )
       if (y_loc(2) .eq. 0.) then
          nm = nm + 1
          rmid(nm) = y_loc(2)
       else if (sign_old .ne. sign_new ) then
          nm = nm + 1
          sign_old = sign_new
          rmid(nm) = y_loc(1)*(x_loc(1) - x_loc(2))/(y_loc(2) - y_loc(1)) + x_loc(1)
       endif
       x_loc(1) = x_loc(2)  ! make new position the old position for next iteration
       y_loc(1) = y_loc(2)
    enddo
    ! find min and mac values
    rmin = 100.
    rmax = -100
    if (nm .eq. 0) then
       ! no interseection found, set values to nonsensical values
       p%x0_min = 1.e10    ! mina 
       p%x0_max = -1.e10
    else
       p%x0_min = minval(rmid(1:nm))
       p%x0_max = maxval(rmid(1:nm))
       p%intersects_midplane = .true.
    endif
    deallocate(rmid)
    return
  end subroutine get_mid_plane

    
  subroutine init_limiter
    ! This routine reads in the data which describes the non-axisymmetric
    ! limiter outline (as regions demarcated by toroidal angle).  This
    ! data is then used by the subroutine INTLIM which tests to see whether
    ! the particle has hit the limiter.
    
    implicit none
    
    real(kind = 8) :: loc_phi
    !     CDUM for skipping comment lines in input file
    character(len = 4) :: cdum
    integer(kind = 4) ::i, j, k, itor, jpol, nmid, npoly, ioerr
    
    nmid = 0 ! for counting number of mid-plane sections for drawing

    if (limiter_directory .eq. '') then
       print *, 'limiter directory no set!'
       limiter_initialized = .false.
       return
    endif
    if (limiter_filename  .eq. '') then
       print *, 'limiter filename no set!'
       limiter_initialized = .false.
       return
    endif
    
    ! read limiter data file
    open (unit=nlimf, file=TRIM(ADJUSTL(limiter_directory))//'/'//limiter_filename, status='old', iostat = ioerr)
    if (ioerr .ne. 0) then
       print *, 'problem opening limiter file : ', TRIM(ADJUSTL(limiter_directory))//'/'//limiter_filename, ' status = ', ioerr
       limiter_initialized = .false.
       return
    endif

    READ(nlimf, *)    ! skip comment line
    READ(nlimf, *) ntorreg
    print *, 'Number of toroidal regions =', ntorreg
    
    if ((ntorreg .gt. ntorrmx) .or. (ntorreg .le. 0)) then
       print *, ' Number of toroidal regions for limiter, ', ntorreg, ' must be between 0 and ', ntorrmx, ', but is not.'
       stop
    endif
    
    !     For number of toroidal regions specified, read starting toroidal angle, 
    !     number of polygons and number of points in each polygon,
    !     then read the actual coordinates.
    
    do itor= 1, ntorreg
       READ(nlimf, '(a4)') cdum
       READ(nlimf, *) loc_phi
       if ((loc_phi .lt. 0.0) .or. (loc_phi .ge. 360.0)) then
          print *, ' Angle defining range of limiter is less than 0 or greater than 360 degrees: ', loc_phi,&
               ' (region ', itor,')'
          stop
       endif
       !       Take input angles as degrees, but convert them now to radians.
       t_regions(itor)%phi_start = loc_phi * dtr
       phi_start(itor) = loc_phi * dtr
       if (itor .eq. 1) then
          if (loc_phi .ne. 0.0) then
             print *, ' Starting toroidal angle of first toroidal region must be 0.0, but is not: ', loc_phi
             stop
          endif
       else
          if ( t_regions(itor)%phi_start .le.  t_regions(itor-1)%phi_start) then
             print *, ' Starting angles of toroidal regions must be in ascending order, but are not for region ', i
             stop
          endif
       endif
       
       ! read the polygons for region itor
       READ(nlimf, *) npoly
       print *, 'toroidal region # ', itor
       print *, 'Phistt = ', loc_phi
       print *, 'Number of polygons =', npoly
       if (npoly .gt. npolymx) then
          print *, ' Number of polygons ', npoly,' exceeds limit of ', npolymx, ' for region ', itor
          stop
       endif
       ! save number of of polygons for region itor
       t_regions(itor)%np = npoly
       ! get number of vertices for each polygon in itor
       READ(nlimf, *) (t_regions(itor)%poly(jpol)%nv, jpol= 1, npoly)
       ! check the range of number of vertices for each polygon
       do jpol=1, npoly
          print *, 'polygon # ', jpol
          if (t_regions(itor)%poly(jpol)%nv .gt. npolyvmx) then
             print *, ' Number of points in polygon, ', t_regions(itor)%poly(jpol)%nv ,&
                  ' exceeds limit of ', npolyvmx, ' for region ', itor, ' polygon ', jpol
             stop
          else
             ! read the vertex data
             print *, 'Polygon for region ', itor, ' polygon ', jpol
             do k=1, t_regions(itor)%poly(jpol)%nv
                ! read the vertex k for polygon jpol, for toroidal region itor
                READ(nlimf, *)  t_regions(itor)%poly(jpol)%x(k),   t_regions(itor)%poly(jpol)%y(k)
             end do
             call print_poly('init_limiter vertex data ', t_regions(itor)%poly(jpol) )
             READ(nlimf, *) ! read blank line separating the polygon data
          endif
       end do
       ! find rmin and rmax at mid-plane for later drawing for the first (main) polygon
       call get_mid_plane(t_regions(itor)%poly(1))
       if ( t_regions(itor)%poly(1)%intersects_midplane ) nmid = nmid + 1
    enddo
    
    t_regions(ntorreg+1)%phi_start = twopi
    t_regions(ntorreg+1)%np  = 0
    phi_start(ntorreg+1) = twopi
    
    close (unit=nlimf)
    limiter_initialized = .true.
    
    ! write the data for plotting
    open(52, file = 'limiter_drawing.data')

    ! R Z poins for the limiter regions 
    do itor= 1, ntorreg
       write (52, *) '-region ', itor, phi_start(i)
       write (52, *) '-region ', itor, t_regions(itor)%phi_start
       
       do k = 1, t_regions(itor)%np
          write (52,*) '--polygon, ', t_regions(itor)%poly(k)%nv  
          do j = 1, t_regions(itor)%poly(k)%nv
             write(52, *) t_regions(itor)%poly(k)%x(j), t_regions(itor)%poly(k)%y(j)
          enddo
       enddo
    enddo

    ! midplane drawing
    write (52,*) '--midplane, ', nmid+1
    write(52, *) '---inner(r,phi)'
    do itor = 1, ntorreg
       if  (t_regions(itor)%poly(1)%intersects_midplane)  then
          write(52, *) t_regions(itor)%poly(1)%x0_min, t_regions(itor)%phi_start
       endif
    enddo
    write(52, *) t_regions(ntorreg)%poly(1)%x0_min, twopi ! to close the drawing
    write(52, *) '---outer(r,phi)'
    do itor = 1, ntorreg
       if  (t_regions(itor)%poly(1)%intersects_midplane) then
          write(52, *) t_regions(itor)%poly(1)%x0_max, t_regions(itor)%phi_start
       endif
    enddo
    write(52, *) t_regions(ntorreg)%poly(1)%x0_max, twopi ! to close the drawing
    close(52)
    return
  end subroutine init_limiter


  
  function  limiter_hit (r, z, phi_v) result(hitlim)
    
    ! This routine checks to see whether a particle has struck the
    ! limiter on its most recent step in its orbit. Return the 
    ! logical variable HITLIM as TRUE if particle has hit the 
    ! limiter, FALSE otherwise.  The particle is considered to have
    ! hit the limiter if it is on the line forming the outline of
    ! the limiter, or is on the side of that line away from where
    ! the plasma is supposed to reside.

    implicit none

    real(kind = 8), intent(in) :: r, z, phi_v
    logical :: hitlim
    integer(kind = 4) :: pos

    ! local variables
    real(kind = 8) :: phi
    integer(kind = 4) :: itor, k, itorreg
    logical :: in_poly

    if (.not. limiter_initialized) then
       print *, 'limiter module not completely initialized '
       hitlim = .true.
       return
    endif
    
    ! set HITLIM to .TRUE. if particle has hit limiter.
    hitlim = .false.
    ! First, compute phi (toroidal position) of particle, modulo 2*pi.
    phi = mod(phi_v, twopi)

    ! Now find the toroidal region in which particle lies
    itor = locate (phi, phi_start, ntorreg+1)
    if ((itor .le. 0) .or. (itor .gt. ntorreg)) then
       print *, ' Particle toroidal angle (mod 2pi) is not in range:', phi, ' r=', r, ' itor=', itor
       hitlim = .true.
       return
    endif
    in_poly =  in_polygon(r, z, t_regions(itor)%poly(1) )
    if ( in_poly ) then
       ! is inside main polygon check the others
       do k=2, t_regions(itor)%np   ! loop over polygons in this regions
          in_poly =  in_polygon(r, z, t_regions(itor)%poly(k) )
          if ( in_poly ) then
             ! inside an interior poly gone -> hit limiter
             hitlim = .true.
             return
          endif
       enddo
       ! did not hit inside polygon
       hitlim = .false.
       return
    else
       ! not inside main polygon
       hitlim = .true.
    endif
    return
    
  end function limiter_hit


  function hit_lim(r) result(hit)
    real(kind = 8), dimension(3), intent(in) :: r
    logical :: hit
    real(kind = 8) :: rt, zt, phit
    
    ! convert to toroidal coordinates
    rt = sqrt(r(1)**2 + r(2)**2)
    zt = r(3)
    phit = pol_angle(r(1:2))

    hit = limiter_hit(rt, zt, phit)
    if (hit .and. print_hit ) then
       print *, '----------------------------------------------------------------------'
       print *, '--- limiter hit at r,z,phi = ', rt, zt, phit/dtr
       print *, '----------------------------------------------------------------------'
    endif
    return
  end function hit_lim
  
  function hit_lim_tor(r) result(hit)
    ! assuming r: r(1) = r, r(2) = phi, r(3) = z
    real(kind = 8), dimension(3), intent(in) :: r
    logical :: hit

    hit = limiter_hit(r(1), r(3), r(2))
    if (hit .and. print_hit ) then
       print *, '----------------------------------------------------------------------'
       print *, '--- limiter hit at r,z,phi = ', r(1), r(3), r(2)/dtr
       print *, '----------------------------------------------------------------------'
    endif
    return
  end function hit_lim_tor
  
end module limiter_mod

