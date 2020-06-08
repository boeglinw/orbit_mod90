	  
	  
	  


subroutine spli2d (tau, gtau, t, n, k, m, work, q, bcoef, iflag)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!  calls bsplvb, banfac/slv
!  this is an extended version of  splint , for the use in tensor prod-
!  uct interpolation.
!
!   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of
!   order  k  with knots  t (i), i=1,..., n + k , which takes on the
!   value  gtau (i,j)  at  tau (i), i=1,..., n ; j=1,..., m .
!
!******  i n p u t  ******
!  tau   array of length  n , containing data point abscissae.
!  a s s u m p t i o n . . .  tau  is strictly increasing
!  gtau(.,j)  corresponding array of length  n , containing data point
!        ordinates, j=1,...,m
!  t     knot sequence, of length  n+k
!  n     number of data points and dimension of spline space  s(k,t)
!  k     order of spline
!  m     number of data sets
!
!******  w o r k   a r e a  ******
!  work  a vector of length  n
!
!******  o u t p u t  ******
!  q     array of order  (n,2*k-1), containing the triangular factoriz-
!        ation of the coefficient matrix of the linear system for the b-
!        coefficients of the spline interpolant.
!           the b-coeffs for the interpolant of an additional data set
!        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
!        be obtained without going through all the calculations in this
!        routine, simply by loading  htau  into  bcoef  and then execut-
!        ing the    call banslv ( q, n, n, 2*k-1, k, bcoef )
!  bcoef the b-coefficients of the interpolant, of length  n
!  iflag an integer indicating success (= 1)  or failure (= 2)
!        the linear system to be solved is (theoretically) invertible if
!        and only if
!              t(i) .lt. tau(i) .lt. tau(i+k),    all i.
!        violation of this condition is certain to lead to  iflag = 2 .
!
!******  m e t h o d  ******
!     the i-th equation of the linear system  a*bcoef = b  for the b-co-
!  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
!  bands (if it is invertible).
!     the matrix  a  is generated row by row and stored, diagonal by di-
!  agonal, in the  c o l u m n s  of the array  q , with the main diag-
!  onal going into column  k .  see comments in the program below.
!     the banded system is then solved by a call to  banfac (which con-
!  structs the triangular factorization for  a  and stores it again in
!   q ), followed by a call to  banslv (which then obtains the solution
!   bcoef  by substitution).
!     banfac  does no pivoting, since the total positivity of the matrix
!  a  makes this unnecessary.
!
!      integer iflag,k,m,n,i,ilp1mx,j,jj,kpkm1,left,np1
!      real bcoef(m,n),gtau(n,m),q(n,7),t(n+k),tau(n),work(n),taui


implicit real (8)(a - h, o - z)  
save  
dimension bcoef (m, n), gtau (n, m), q (n, 2 * k - 1), t (n + k), &
 tau (n), work (n)
!
nnn = 1  
np1 = n + 1  
kpkm1 = 2 * k - 1  
left = k  
!
!  ***   loop over  i  to construct the  n  interpolation equations
do 30 i = 1, n  
   iindex = i  
   taui = tau (iindex)  
   ilp1mx = min (iindex + k, np1)  
!        *** zero out all entries in row  i  of  a (in the 2k-1 bands)
   do 13 j = 1, kpkm1  
   13    q (iindex, j) = 0.  
!        *** find  left  in the closed interval (i,i+k-1) such that
!                t(left) .le. tau(i) .lt. t(left+1)
!        matrix is singular if this is not possible
   left = max (left, i)  
   if (taui.lt.t (left) ) goto 998  
   15    if (taui.lt.t (left + 1) ) goto 16  
   left = left + 1  
   if (left.lt.ilp1mx) goto 15  
   left = left - 1  
   if (taui.gt.t (left + 1) ) goto 998  
!        *** the i-th equation enforces interpolation at taui, hence
!        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
!        left-k+1,...,left actually might be nonzero. these  k  numbers
!        are returned, in  work  (used for temp.storage here), by the
!        following
   16    call bsplvb (t, k, nnn, taui, left, work)  
!        we therefore want  work(j) = b(left-k+j)(taui) to go into
!        a(i,left-k+j), i.e., into  q(i,left-i+j), since the i-th row of
!        a  is so stored in the i-th row of  q  that the (i,i)-entry of
!        a  goes into the  k-th  entry of  q.
   jj = left - iindex  
   do 29 j = 1, k  
      jj = jj + 1  
      q (iindex, jj) = work (j)  
   29    end do  
   30 end do  
!
!     ***obtain factorization of  a  , stored again in  q.
call banfac (q, n, n, kpkm1, k, iflag)  
goto (40, 999), iflag  
!     *** solve  a*bcoef = gtau  by backsubstitution
   40 do 50 j = 1, m  
   do 41 i = 1, n  
   41    work (i) = gtau (i, j)  
   call banslv (q, n, n, kpkm1, k, work)  
   do 50 i = 1, n  
   50 bcoef (j, i) = work (i)  
return  
  998 iflag = 2  
  999 print 699  
  699 format(41h linear system in  splint  not invertible)  
return  
end subroutine spli2d
