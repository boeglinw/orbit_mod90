
!
!   This routine is required if the CVS revision numbers are to
!   survive an optimization.
!
!
!   $Date: 1997/04/05 01:43:16 $ $Author: peng $
!



subroutine spline_rev (i)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
CHARACTER (len=100) :: opt  
character (len=10) :: s  
if (i.eq.0) s = '@(#)$RCSfile: spline.f,v $ $Revision: 2.1 $\000'  
return  
end subroutine spline_rev
