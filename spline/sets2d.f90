



subroutine sets2d (s, cs, x, nx, bkx, lx, y, ny, bky, ly, wk, ier)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!-----------------------------------------------------------------------
!--  S.Thompson  92/05/18
!--    Bicubic spline routines.
!--    Put together with routines from E.Solano.
!--  SMWolfe     93/12/17
!--    Modifed to avoid over-writing the original array.
!--              94/04/28
!--    Updated.
!-----------------------------------------------------------------------
!  Inputs:
!
!      s     - nx by ny array containing the function values at (x,y).
!              This is a 1-d array, k=k=(i-1)*ny+j.
!
!      x, y  - (x,y) location, arrays of length nx and ny.
!
!  Outputs:
!
!      cs    - array of spline coefficients of dimension (kubicx,
!              lubicx,kubicy,lubicy).
!
!      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1.
!
!      lx, ly -   number of terms in bkx and bky.
!
!      ier   - rror parameter.
!
!  Work arrays:
!
!      wk    - of dimension at least nx by ny.
!-----------------------------------------------------------------------

use orbit_parameters

implicit real (8)(a - h, o - z)  
save  



!
!      dimension s(1), x(nx), y(ny), wk(nx,ny),
!     .          xknot(kubicx + nw), yknot(kubicy + nh),
!     .          cs(kubicx, lubicx, kubicy, lubicy),
!     .          bkx(lubicx + 1), bky(lubicy + 1)
dimension s (1), x (nx), y (ny), wk (nx, ny), cs (kubicx, nx - &
 kubicx + 1, kubicy, ny - kubicy + 1), bkx (nx - kubicx + 2), &
 bky (ny - kubicy + 2), xknot (nxknot), yknot (nyknot)
!     .          xknot(kubicx + nx), yknot(kubicy + ny),
!
!  Set up knots:
!
call eknot (nx, x, kubicx, xknot) 		  
call eknot (ny, y, kubicy, yknot) 			  
!
!  Save the original, use the work array
!
do 10 i = 1, nx  
   do 10 j = 1, ny  
      k = (i - 1) * ny + j  
   10 wk (i, j) = s (k)  
!
!  Calculate spline coefficients:
!
call spl2bc (x, y, nx, ny, xknot, yknot, wk) 	  
!
!  Coefficients stored in bkx, bky, and c:
!
call spl2pp (nx, ny, xknot, yknot, wk, bkx, lx, bky, ly, cs)  
!
return  
end subroutine sets2d
