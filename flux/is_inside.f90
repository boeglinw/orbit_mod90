

function is_inside (x, y, nlim, xlim, ylim)  
  !**********************************************************************
  !                                                               
  !                                                               
  !     SUBPROGRAM DESCRIPTION:                                   
  !          is_inside determines whether point (x,y) is  
  !          inside or outside of the boundary given by the nlim
  !          number of points xlim, ylim                          
  !                                                                  
  !                                                                  
  !     CALLING ARGUMENTS:                                           
  !       x...............x-coordinate of point (input)                            
  !       y...............y-coordinate of point (input)                            
  !       nlim............number of boundary points (input)           
  !       xlim............x coordinates of boundary (input)           
  !       ylim............y coordinates of boundar  (input)           
  !                                                                  
  !     REFERENCES:                                                  
  !          (1)  uses ray casting
  !                                                                  
  !     RECORD OF MODIFICATION:                                      
  !          created 06/03/2020... W. Boeglin                  
  !                                                                  
  !                                                                  
  !                                                                  
  !**********************************************************************
  
  implicit none

  real(kind = 8) :: x
  real(kind = 8) :: y   
  integer(kind = 4) :: nlim 
  real(kind = 8), dimension(nlim) :: xlim
  real(kind = 8), dimension(nlim) :: ylim
  
  real(kind = 8) :: s, t, di, f   

  logical :: is_inside

  integer(kind = 4) :: i, j, k, ncross
  
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
  is_inside = mod(ncross,2) .eq. 1 ! inside if odd number of crossings for cast ray
  return  
end function is_inside

