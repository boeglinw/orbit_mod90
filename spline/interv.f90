



subroutine interv (xt, lxt, x, left, mflag)  
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
!omputes  left = max( i ; 1 .le. i .le. lxt  .and.  xt(i) .le. x )  .
!
!******  i n p u t  ******
!  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
!  lxt.....number of terms in the sequence  xt .
!  x.....the point whose location with respect to the sequence  xt  is
!        to be determined.
!
!******  o u t p u t  ******
!  left, mflag.....both integers, whose value is
!
!   1     -1      if               x .lt.  xt(1)
!   i      0      if   xt(i)  .le. x .lt. xt(i+1)
!  lxt     1      if  xt(lxt) .le. x
!
!        in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
!        indicates that  x  lies outside the halfopen interval
!        xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
!        interval is due to the decision to make all pp functions cont-
!        inuous from the right.
!
!******  m e t h o d  ******
!  the program is designed to be efficient in the common situation that
!  it is called repeatedly, with  x  taken from an increasing or decrea-
!  sing sequence. this will happen, e.g., when a pp function is to be
!  graphed. the first guess for  left  is therefore taken to be the val-
!  ue returned at the previous call and stored in the  l o c a l  varia-
!  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
!  essary since the present call may have nothing to do with the previ-
!  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!  ilo  and are done after just three comparisons.
!     otherwise, we repeatedly double the difference  istep = ihi - ilo
!  while also moving  ilo  and  ihi  in the direction of  x , until
!                      xt(ilo) .le. x .lt. xt(ihi) ,
!  after which we use bisection to get, in addition, ilo+1 = ihi .
!  left = ilo  is then returned.
!

implicit real (8)(a - h, o - z)  
save  
integer :: left, lxt, mflag, ihi, ilo, istep, middle  
real (8) :: x  
dimension xt (lxt)  
data ilo / 1 /  
!     save ilo  (a valid fortran statement in the new 1977 standard)
ihi = ilo + 1  
if (ihi.lt.lxt) goto 20  
if (x.ge.xt (lxt) ) goto 110  
if (lxt.le.1) goto 90  
ilo = lxt - 1  
ihi = lxt  
!
   20 if (x.ge.xt (ihi) ) goto 40  
if (x.ge.xt (ilo) ) goto 100  
!
!              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
   30 istep = 1  
   31 ihi = ilo  
ilo = ihi - istep  
if (ilo.le.1) goto 35  
if (x.ge.xt (ilo) ) goto 50  
istep = istep * 2  
goto 31  
   35 ilo = 1  
if (x.lt.xt (1) ) goto 90  
goto 50  
!              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1  
   41 ilo = ihi  
ihi = ilo + istep  
if (ihi.ge.lxt) goto 45  
if (x.lt.xt (ihi) ) goto 50  
istep = istep * 2  
goto 41  
   45 if (x.ge.xt (lxt) ) goto 110  
ihi = lxt  
!
!           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi) / 2  
if (middle.eq.ilo) goto 100  
!     note. it is assumed that middle = ilo in case ihi = ilo+1 .
if (x.lt.xt (middle) ) goto 53  
ilo = middle  
goto 50  
   53 ihi = middle  
goto 50  
!**** set output and return.
   90 mflag = - 1  
left = 1  
return  
  100 mflag = 0  
left = ilo  
return  
  110 mflag = 1  
left = lxt  
return  
end subroutine interv
