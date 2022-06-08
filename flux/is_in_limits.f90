


subroutine is_in_limits(np, xp, yp, nlim, xlim, ylim, inside)
  !**********************************************************************
  !                                                               
  !                                                               
  !     SUBPROGRAM DESCRIPTION:                                   
  !          is_inside determines whether array of points (xp,yp) are  
  !          inside or outside of the boundary given by the nlim
  !          number of points xlim, ylim                          
  !                                                                  
  !                                                                  
  !     CALLING ARGUMENTS:
  !       np...............number of points                                        
  !       xp...............x-coordinate of points (input)                            
  !       yp...............y-coordinate of points (input)                            
  !       nlim............number of boundary points (input)           
  !       xlim............x coordinates of boundary (input)           
  !       ylim............y coordinates of boundar  (input)           
  !                                                                  
  !     REFERENCES:                                                  
  !          (1)  uses ray casting
  !                                                                  
  !     RECORD OF MODIFICATION:                                      
  !          created 06/08/2022... W. Boeglin                  
  !                                                                  
  !                                                                  
  !                                                                  
  !**********************************************************************
  
  implicit none

  ! data points
  integer(kind = 4) :: np 
  real(kind = 8), dimension(np), intent(in) :: xp
  real(kind = 8), dimension(np), intent(in) :: yp
  ! results  
  logical, dimension(np), intent(out)::inside 
  
  ! boundary points
  integer(kind = 4) :: nlim 
  real(kind = 8), dimension(nlim), intent(in) :: xlim
  real(kind = 8), dimension(nlim), intent(in) :: ylim
  
  real(kind = 8) :: x, y
  
  real(kind = 8) :: s, t, di, f   
  
  integer :: i

  integer(kind = 4) ::  k, ncross

  do i = 1, np
      x = xp(i)
      y = yp(i)
      ! cast the horizontal ray and count the number of crossings, if odd then the point is inside
      ncross = 0  
      do k = 1, nlim - 1
         if ( (ylim (k) .lt.y ) .and. (ylim (k + 1) .lt.y ) ) cycle
         if (x .eq. xlim (k) ) cycle  
         t = x - xlim (k)  
         s = xlim (k + 1) - x  
         if ( (t * s) .lt.0.) cycle  
         di = (ylim (k + 1) - ylim (k) ) / (xlim (k + 1) - xlim (k) )
         f = ylim (k) + di * (x - xlim (k) )  
         if (f.lt.y ) cycle  
         ncross = ncross + 1
      end do
      inside(i) = (mod(ncross,2) .eq. 1) ! inside if odd number of crossings for cast ray
  enddo
  return  
end subroutine is_in_limits

