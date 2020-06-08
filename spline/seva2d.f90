
subroutine seva2d (bkx, lx, bky, ly, cs, nx, ny, xl, yl, fs, ier, &
 icalc)
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
!      cs       - array of spline coefficients of dimension (kubicx,
!                 lubicx,kubicy,lubicy) from sets2d.
!
!      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1
!                 sets2d.
!
!      lx, ly   - number of terms in bkx and bky from sets2d.
!
!      xl, yl   - the point at which interpolations are desired.
!
!      nx, ny   - grid dimensions
!
!  Outputs:
!
!      fs       - vector containing results depending on icalc:
!                 icalc              fs
!                   1                f
!                   2                fx
!                   3                fy
!                   4                fxy
!                   5                fxx
!                   6                fyy
!
!      ier      - error parameter.
!
!-----------------------------------------------------------------------

  use orbit_parameters

implicit real (8)(a - h, o - z)  
save  


!
integer :: ier, lx, ly  
!      real cs(kubicx,lubicx,kubicy,lubicy),xl,yl,fs(6),bkx(1),bky(1)
real (8) :: xl, yl, fs (6), bkx (1), bky (1)  
dimension cs (kubicx, nx - kubicx + 1, kubicy, ny - kubicy + 1)  
!
!  Local Variable Specifications:
!
dimension work0 (4), work1 (4), work2 (4)  
data n00 / 0 /, n11 / 1 /, n22 / 2 /  
!	  integer n00, n11, n22
!
!  Evaluate function and its partial derivatives at (XL, YL):
!
!  First do all the lookup and interpolation stuff.
!  This is the most time consuming part of the evaluation, so
!  don't do more than needed.
!
call interv (bky, ly, yl, lef, mflag)  
call interv (bkx, lx, xl, ibk, ndummy)  
h = xl - bkx (ibk)  
do 41 jj = 1, 4  
   work0 (jj) = ppvalw (cs (1, ibk, jj, lef), h, n00)  
   if (icalc.eq.1) goto 41  
   work1 (jj) = ppvalw (cs (1, ibk, jj, lef), h, n11)  
   if (icalc.le.4) goto 41  
   work2 (jj) = ppvalw (cs (1, ibk, jj, lef), h, n22)  
   41 end do  
h = yl - bky (lef)  
fs (1) = ppvalw (work0, h, n00)  
if (icalc.eq.1) return  
fs (2) = ppvalw (work1, h, n00)  
if (icalc.eq.2) return  
fs (3) = ppvalw (work0, h, n11)  
if (icalc.eq.3) return  
fs (4) = ppvalw (work1, h, n11)  
if (icalc.eq.4) return  
fs (5) = ppvalw (work2, h, n00)  
if (icalc.eq.5) return  
fs (6) = ppvalw (work0, h, n22)  
!
return  
end subroutine seva2d
