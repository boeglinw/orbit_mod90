



subroutine bsplvb (t, jhigh, index, x, left, biatx)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!  calculates the value of all possibly nonzero b-splines at  x  of orde
!
!               jout  =  max( jhigh , (j+1)*(index-1) )
!
!  with knot sequence  t .
!
!******  i n p u t  ******
!  t.....knot sequence, of length  left + jout  , assumed to be nonde-
!        creasing.  a s s u m p t i o n . . . .
!                       t(left)  .lt.  t(left + 1)   .
!   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
!  jhigh,
!  index.....integers which determine the order  jout = max(jhigh,
!        (j+1)*(index-1))  of the b-splines whose values at  x  are to
!        be returned.  index  is used to avoid recalculations when seve-
!        ral columns of the triangular array of b-spline values are nee-
!        ded (e.g., in  bvalue  or in  bsplvd ). precisely,
!                     if  index = 1 ,
!        the calculation starts from scratch and the entire triangular
!        array of b-spline values of orders 1,2,...,jhigh  is generated
!        order by order , i.e., column by column .
!                     if  index = 2 ,
!        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
!        nerated, the assumption being that  biatx , j , deltal , deltar
!        are, on entry, as they were on exit at the previous call.
!           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
!        the next column of b-spline values is generated.
!
!  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
!        posed arbitrarily by the dimension statement for  deltal  and
!        deltar  below, but is  n o w h e r e  c h e c k e d  for .
!
!  x.....the point at which the b-splines are to be evaluated.
!  left.....an integer chosen (usually) so that
!                  t(left) .le. x .le. t(left+1)  .
!
!******  o u t p u t  ******
!  biatx.....array of length  jout , with  biatx(i)  containing the val-
!        ue at  x  of the polynomial of order  jout  which agrees with
!        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
!        t(left+1)) .
!
!******  m e t h o d  ******
!  the recurrence relation
!
!                       x - t(i)              t(i+j+1) - x
!     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
!                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
!
!  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
!  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
!  b(left,j)(x), storing the new values in  biatx  over the old. the
!  facts that
!            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
!  and that
!            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
!  are used. the particular organization of the calculations follows al-
!  gorithm  (8)  in chapter x of the text.
!

implicit real (8)(a - h, o - z)  
save  
parameter (jmax = 4)  
integer :: index, jhigh, left, i, j, jp1  
real (8) :: x, saved, term  
!      real biatx(jhigh),t(1),x,
dimension deltal (jmax), deltar (jmax)  
dimension biatx (jhigh), t (left + jhigh)  
!  current fortran standard makes it impossible to specify the length of
!  t and of biatx  precisely without the introduction of otherwise
!  superfluous additional arguments.
data j / 1 /  
!      save j,deltal,deltar  ! (valid in fortran 77)
!
goto (10, 20), index  
   10 j = 1  
biatx (1) = 1.  
if (j.ge.jhigh) goto 99  
!
   20 jp1 = j + 1  
deltar (j) = t (left + j) - x  
deltal (j) = x - t (left + 1 - j)  
saved = 0.  
do 26 i = 1, j  
   term = biatx (i) / (deltar (i) + deltal (jp1 - i) )  
   biatx (i) = saved+deltar (i) * term  
   26 saved = deltal (jp1 - i) * term  
biatx (jp1) = saved  
j = jp1  
if (j.lt.jhigh) goto 20  
!
   99 return  
end subroutine bsplvb
