




subroutine spl2bc (rgrid, zgrid, mw, mh, rknot, zknot, copy)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
! calculates the b-spline coeficients

use orbit_parameters_mod

implicit real (8)(a - h, o - z)  
save  


! changed WB 2018, to using standard size

! parameter (nw = 129, nh = 129, krord = 4, kzord = 4)  

parameter (krord = 4, kzord = 4)  

dimension rgrid (mw), zgrid (mh)  

!      dimension rknot(nw+krord),zknot(nh+kzord),copy(mw,mh)


dimension rknot (mw + krord), zknot (mh + kzord), copy (mw, mh)  
!------------------------------------------------------------------
!-- change dimension of work2 and work3 from nw to nh            --
!-- to ensure the cases when nh > nw     ll, 93/04/01            --
!------------------------------------------------------------------
!      dimension work1(mw,mh),work2(mh),work3(mh,2*krord-1)

dimension work1 (nw, nh), work2 (nh), work3 (nh, 2 * krord-1)  
call spli2d (rgrid, copy, rknot, mw, krord, mh, work2, work3, &
 work1, iflag)

if (iflag.ne.1) print * , ' error in first spli2d, iflag=', iflag  
call spli2d (zgrid, work1, zknot, mh, kzord, mw, work2, work3, &
 copy, iflag)
if (iflag.ne.1) print * , ' error in second spli2d, iflag=', &
 iflag
return  
end subroutine spl2bc
