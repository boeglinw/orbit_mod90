




subroutine bspp2d (t, bcoef, n, k, m, scrtch, break, coef, l)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!  calls  bsplvb
!  this is an extended version of  bsplpp  for use with tensor products
!
!onverts the b-representation  t, bcoef(.,j), n, k  of some spline into
!  its pp-representation  break, coef(j,.,.), l, k ; j=1, ..., m  .
!
!******  i n p u t  ******
!  t     knot sequence, of length  n+k
!  bcoef(.,j) b-spline coefficient sequence, of length  n ;j=1,...,m
!  n     length of  bcoef  and  dimension of spline space  s(k,t)
!  k     order of the spline
!  m     number of data sets
!
!******  w o r k   a r e a  ******
!  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of
!        the spline and its  k-1  derivatives   for each of the m sets
!
!******  o u t p u t  ******
!  break breakpoint sequence, of length  l+1, contains (in increasing
!        order) the distinct points in the sequence  t(k), ..., t(n+1)
!  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der-
!        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m
!  l     number of polynomial pieces which make up the spline in the
!        interval  (t(k), t(n+1))
!
!******  m e t h o d  ******
!     for each breakpoint interval, the  k  relevant b-coeffs of the
!  spline are found and then differenced repeatedly to get the b-coeffs
!  of all the derivatives of the spline on that interval. the spline and
!  its first  k-1  derivatives are then evaluated at the left end
!  point of that interval, using  bsplvb  repeatedly to obtain the val-
!  ues of all b-splines of the appropriate order at that point.
!

implicit real (8)(a - h, o - z)  
save  
parameter (kmax = 4)  
integer :: k, l, m, n, i, j, jp1, kmj, left  
dimension bcoef (n, m), break (1), coef (m, k, 1), scrtch (k, k, &
 m), t (1), biatx (kmax)
real (8) :: diff, fkmj, sum  
!
n11 = 1  
n22 = 2  
l = 0  
break (1) = t (k)  
do 50 left = k, n  
!        find the next nontrivial knot interval.
   if (t (left + 1) .eq.t (left) ) goto 50  
   l = l + 1  
   break (l + 1) = t (left + 1)  
   if (k.gt.1) goto 9  
   do 5 mm = 1, m  
    5    coef (mm, 1, l) = bcoef (left, mm)  
   goto 50  
!        store the k b-spline coeff.s relevant to current knot interval
!        in  scrtch(.,1) .
    9    do 10 i = 1, k  
      do 10 mm = 1, m  
   10    scrtch (i, 1, mm) = bcoef (left - k + i, mm)  
!        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
!        current knot interval for the j-th derivative by differencing
!        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
   do 20 jp1 = 2, k  
      j = jp1 - 1  
      kmj = k - j  
      fkmj = float (kmj)  
      do 20 i = 1, kmj  
         diff = (t (left + i) - t (left + i - kmj) ) / fkmj  
         if (diff.le.0.) goto 20  
         do 15 mm = 1, m  
   15          scrtch (i, jp1, mm) = (scrtch (i + 1, j, mm) - scrtch (i, &
          j, mm) ) / diff
   20    continue  
!        starting with the one b-spline of order 1 not zero at t(left),
!        find the values at t(left) of the j+1 b-splines of order j+1
!        not identically zero there from those of order j, then combine
!        with the b-spline coeff.s found earlier to compute the (k-j)-
!        th derivative at t(left) of the given spline.
   call bsplvb (t, n11, n11, t (left), left, biatx)  
   do 25 mm = 1, m  
   25    coef (mm, k, l) = scrtch (1, k, mm)  
   do 30 jp1 = 2, k  
      call bsplvb (t, jp1, n22, t (left), left, biatx)  
      kmj = k + 1 - jp1  
      do 30 mm = 1, m  
         sum = 0.  
         do 28 i = 1, jp1  
   28          sum = biatx (i) * scrtch (i, kmj, mm) + sum  
   30    coef (mm, kmj, l) = sum  
   50 end do  
return  
end subroutine bspp2d
