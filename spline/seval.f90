



real(8) function seval (n, u, x, y, b, c, d)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************

implicit real (8)(a - h, o - z)  
save  
integer :: n  
real (8) :: u, x (n), y (n), b (n), c (n), d (n)  
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by zpline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
integer :: i, j, k  
real (8) :: dx  
data i / 1 /  
if (i.ge.n) i = 1  
if (u.lt.x (i) ) goto 10  
if (u.le.x (i + 1) ) goto 30  
!
!  binary search
!
   10 i = 1  
j = n + 1  
   20 k = (i + j) / 2  
if (u.lt.x (k) ) j = k  
if (u.ge.x (k) ) i = k  
if (j.gt.i + 1) goto 20  
!
!  evaluate spline
!
   30 dx = u - x (i)  
seval = y (i) + dx * (b (i) + dx * (c (i) + dx * d (i) ) )  
return  
end function seval
