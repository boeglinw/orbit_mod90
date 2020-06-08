



function ppvalw (coef, x, jd)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!-----------------------------------------------------------------------
!  Modified for optimization by S.J. Thompson, 30-Aug-1993
!  Revised to eliminate call to interv by S.M.Wolfe, 17-Dec-1993
!          and to use ASF's for evaluation
!  This routine performs only the innermost guts of the spline evaluatio
!  Assumes k=4 (cubic spline only). No other cases considered.
! does not call  interv
!alculates value at  x  of  jd-th derivative of pp fct from pp-repr
!
!******  i n p u t  ****** to PPVALU, on which this is based.
!  break, coef, l, k.....forms the pp-representation of the function  f
!        to be evaluated. specifically, the j-th derivative of  f  is
!        given by
!
!     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
!                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
!
!        with  h = x - break(i),  and
!
!       i  =  max( 1 , max( j ;  break(j) .le. x , 1 .le. j .le. l ) ).
!
!  x.....the point at which to evaluate.
!        as used here, x is the distance from the break, not the absolut
!        position.
!  jd.....integer*4 giving the order of the derivative to be evaluat-
!        ed.  a s s u m e d  to be zero or positive.
!
!******  o u t p u t  ******
!  ppvalw.....the value of the (jd)-th derivative of  f  at  x.
!
!******  m e t h o d  ******
!     the interval index  i , appropriate for  x , is found through a
!  call to  interv . the formula above for the  jd-th derivative
!  of  f  is then evaluated (by nested multipication).
!
!-----------------------------------------------------------------------
!   Variable declarations.
!-----------------------------------------------------------------------


implicit real (8)(a - h, o - z)  
save  
real (8) :: ppvalw, x  
dimension coef (4)  
!----------------------------------------------------------------------
! ASF's may be slightly more efficient than the alternative
!----------------------------------------------------------------------
d2 (xx) = coef (4) * xx + coef (3)  
d1 (xx) = (coef (4) * xx / 2. + coef (3) ) * xx + coef (2)  
d0 (xx) = ( (coef (4) * xx / 3. + coef (3) ) * xx / 2. + coef (2) &
 ) * xx + coef (1)
!-----------------------------------------------------------------------
!   Derivatives of order k or higher are identically zero.
!-----------------------------------------------------------------------
!   Evaluate jd-th derivative of i-th polynomial piece at x .
!-----------------------------------------------------------------------
goto (1, 2, 3) jd+1  
ppvalw = 0.  
print * , 'Error (ppvalw): JD must be 0, 1, or 2.'  
print * , 'Execution terminated.'  
return  
                     ! k = 4 , jd = 0
    1 ppvalw = d0 (x) 	  
return  
                     ! k = 4 , jd = 1
    2 ppvalw = d1 (x) 	  
return  
                     ! k = 4 , jd = 2
    3 ppvalw = d2 (x) 	  
return  
end function ppvalw
