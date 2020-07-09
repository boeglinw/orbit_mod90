

subroutine checklim (nx, ny, x, y, nlim, xlim, ylim, inside )  
  !**********************************************************************
  !                                                               
  !                                                               
  !     SUBPROGRAM DESCRIPTION:                                   
  !          check_lim determines whether points on an (x,y) grid are  
  !          inside or outside of the boundary given by the nlim
  !          number of points xlim, ylim, e.g. the last 
  !          flux surface                                         
  !                                                                  
  !                                                                  
  !     CALLING ARGUMENTS:                                           
  !       nx..............dimension of x (input)                     
  !       ny..............dimension of y (input)                     
  !       x...............r grid (input)                             
  !       y...............z grid (input)                             
  !       nlim............number of limiter points (input)           
  !       xlim............r coordinates of limiter (input)           
  !       ylim............z coordinates of limiter (input)           
  !
  !       inside............1 if inside and 0 otherwise (output)     
  !                                                                  
  !     REFERENCES:                                                  
  !          (1)                                                     
  !          (2)                                                     
  !                                                                  
  !     RECORD OF MODIFICATION:
  !          based on old zlim PPPL code
  !          06/22/2018...changed to f90 W. Boeglin                  
  !                                                                  
  !                                                                  
  !                                                                  
  !**********************************************************************
  
  implicit none

  integer(kind = 4), dimension(nx*ny) :: inside
  integer(kind = 4) :: nx, ny, nlim 
  real(kind = 8), dimension(nlim) :: xlim
  real(kind = 8), dimension(nlim) :: ylim
  real(kind = 8), dimension(nx) :: x
  real(kind = 8), dimension(ny) :: y   
  
  real(kind = 8) :: s, t, di, f   

  logical :: is_inside
  logical :: first   = .true.

  integer(kind = 4) :: i, j, k, kk, ncross, mcross, mcross1

  kk = 0  
  do i = 1, nx  
     do j = 1, ny
        kk = kk + 1
        if ( is_inside(x(i), y(j), nlim, xlim, ylim) ) then
           inside(kk) = 1
        else
           inside(kk) = 0
        endif
     enddo
  enddo  
  return  
end subroutine checklim

