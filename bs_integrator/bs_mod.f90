module bs_mod
  use constants_and_masses_mod
  use helper_functions_mod
  use limiter_mod
  use em_fields_mod
  use control_mod
  implicit none

  integer(kind = 4), parameter :: nmax=10
  integer(kind = 4), parameter :: imax=11
  integer(kind = 4), parameter :: nuse=7

  real(kind = 8),  parameter :: one = 1.d0
  real(kind = 8),  parameter :: shrink = .95d0
  real(kind = 8),  parameter :: grow = 1.2d0

  ! tolerance for BS step
  real(kind = 8) :: err_tol = 1.0d-5

contains


  function bs_push( rs, vs, dt, nsteps, track) result(n_calc)
    ! calculates trajectory and returns the actual number of steps calculated, this should
    ! be smaller than nsteps if a limiter is hit and check_limiter is set to .true.
    ! interfaces
    ! initial position

    real(kind = 8), dimension(3), intent(in):: rs
    ! intital velocity
    real(kind = 8), dimension(3), intent(in):: vs
    ! time step
    real(kind = 8) :: dt
    ! total number of steps
    integer(kind = 4), intent(in) ::  Nsteps
    integer(kind = 4) :: n_calc ! calculated steps
    ! track data
    real(kind = 8), dimension(nsteps,6), intent(out) :: track
    
    ! local variables
    ! loop variables
    integer(kind = 4) i, j, k, flag

    real(kind = 8) :: t ! total time
    real(kind = 8), dimension(6) :: r

    ! save first position
    r(1:3) = rs
    r(4:6) = vs
    ! store initial position
    track(1,1:6) = r
    do i = 1, Nsteps-1
       ! calculate step
       call bs_ode( t, dt, r, 6, lorentz, flag)
       ! if the flag is set print it
       if (flag .ne. 0) print *, 'bs_pus: bs_ode flag = ', flag, ' result probably bad'
       ! save position and velocity data
       track(i+1,1:6) = r
       if (check_limiter) then   ! if selected check a limiter hit
          if ( hit_lim(r(1:3)) ) then
             n_calc = i + 1          ! adjust the total number of steps calculated
             return              ! if a limiter is hit return
          endif
       endif
    enddo
    n_calc = Nsteps 
  end function bs_push


  function bs_push_tor( rs, vs, dt, nsteps, track) result(n_calc)
    ! bs_push but in toroidal coordinates
    ! this is the original calcualtion method used in orbit. There is no difference between cartesian and toroidal
    ! calculations. It is here just for historical reasons.
    ! calculates trajectory and returns the actual number of steps calculated, this should
    ! be smaller than nsteps if a limiter is hit and check_limiter is set to .true.
    ! interfaces
    ! initial position
    ! rs : r phi z
    ! vs : vr vphi vz

    real(kind = 8), dimension(3), intent(in):: rs
    ! intital velocity
    real(kind = 8), dimension(3), intent(in):: vs
    ! time step
    real(kind = 8) :: dt
    ! total number of steps
    integer(kind = 4), intent(in) ::  Nsteps
    integer(kind = 4) :: n_calc ! calculated steps
    ! track data
    real(kind = 8), dimension(nsteps,6), intent(out) :: track
    
    ! local variables
    ! loop variables
    integer(kind = 4) i, j, k, flag

    real(kind = 8) :: t ! total time
    real(kind = 8), dimension(6) :: r

    ! save first position
    r(1:3) = rs
    r(4:6) = vs
    ! store initial position
    track(1,1:6) = r
    do i = 1, Nsteps -1
       ! calculate step
       call bs_ode( t, dt, r, 6, lorentz_tor, flag)
       ! if the flag is set print it
       if (flag .ne. 0) print *, 'bs_pus: bs_ode flag = ', flag, ' result probably bad'
       ! save position and velocity data
       track(i+1,1:6) = r
       if (check_limiter) then   ! if selected check a limiter hit
          if ( hit_lim_tor(r(1:3)) ) then
             n_calc = i+1          ! adjust the total number of steps calculated
             return              ! if a limiter is hit return
          endif
       endif
    enddo
    n_calc = Nsteps
  end function bs_push_tor
  

  
  subroutine lorentz(t, r, drdt, nn)
    ! generalized position r: x,y,z, vz, vy, vz
    ! this function needs to be provided
    ! interface to the optional b_field routine
    ! interface to the global (user supplied and linked) B-field routine mag_field
    integer(kind = 4), intent(in) :: nn
    real(kind = 8), intent(in) :: t
    real(kind = 8), dimension(nn), intent(in) :: r
    real(kind = 8), dimension(nn), intent(out) :: drdt

    ! local field
    real(kind = 8), dimension(3) :: B 
    B = bfield3(r(1:3))

    ! calculate derivatives from eq. of motion    
    drdt(1) = r(4)
    drdt(2) = r(5)
    drdt(3) = r(6)
    drdt(4:6) = q_over_m * cross_product(r(4:6), B)
    return
  end subroutine lorentz

  
  subroutine lorentz_tor(t,rv,rvpr,n)

    ! eq. of motion in toroidal coords from orbit3
    
    ! generalized position r: x,y,z, vz, vy, vz
    ! this function needs to be provided

    ! interface to the global (user supplied and linked) B-field routine mag_field
    ! bfield returns vector bf = ( bpolr, bpolz, bphi, btotal, psi_rel)
    integer(kind = 4), intent(in) :: n
    real(kind = 8), intent(in) :: t
    real(kind = 8), dimension(n), intent(in) :: rv
    real(kind = 8), dimension(n), intent(out) :: rvpr

    ! local field
    real(kind = 8), dimension(5) :: bo, bf
    real(kind = 8), dimension(4) :: r ! position

    ! assignment
    !     rv(1) : r
    !     rv(2) : phi
    !     rv(3) : z          
    !     b(1)  : br
    !     b(2)  : bphi
    !     b(3)  : bz
    
    bf = bfield(rv(1), rv(3)) ! time reversal is handled in the em_fields_mod
    
    bo(1) = bf(1) * omega
    bo(2) = bf(3) * omega
    bo(3) = bf(2) * omega

    ! time derivatives
    ! remember |v| = 1
    rvpr(1)=rv(4)
    rvpr(2)=rv(5)/rv(1)
    rvpr(3)=rv(6)
    rvpr(4)=rv(5)*rv(5)/rv(1)+rv(5)*bo(3)-rv(6)*bo(2)
    rvpr(5)=-rv(4)*rv(5)/rv(1)+rv(6)*bo(1)-rv(4)*bo(3)
    rvpr(6)=rv(4)*bo(2)-rv(5)*bo(1)
    return

  end subroutine lorentz_tor

! Bulirsch-Stoer Ordinary Differental Equation driver
! with Runge Kutta alternative when BS fails
!
! y is changed upon return
! s is changed upon return
!
! interval check alrgorithm
!
! x1 < x < x2 : (x - x1)*(x - x2) < 0 : inside
!               (x - x1)*(x - x2) = 0 : at one of the limits
!               (x - x1)*(x - x2) > 0 : outside

  subroutine bs_ode(s, ds, y, neqn, derivs, iflag)
     ! s       : independent variable,changed on return
     ! ds      : stepsize
     ! y       :  array of dependent values changed after step s is taken
     ! neqn    :  # of equations to be integrated
     ! derivs  : user supplied derivative subroutine
     ! relerr  :error tolerances
     ! iflag   : return status 0 = ok, 5 = not converged, 6 = tolerance <=0.0

    interface
       subroutine derivs(xx, yy, dyydxx, nn)
         integer(kind = 4), intent(in) :: nn
         real(kind = 8), intent(in) :: xx
         real(kind = 8), dimension(nn), intent(in) :: yy
         real(kind = 8), dimension(nn), intent(out) :: dyydxx
       end subroutine derivs
    end interface

    integer(kind = 4), parameter :: maxstep = 25
    integer(kind = 4), parameter :: max_bsstep_tries = 10
    integer(kind = 4), parameter :: maxrk = 3 

    ! arguments
    integer(kind = 4), intent(in) :: neqn

    real(kind = 8), intent(inout) :: s
    real(kind = 8), intent(in)  :: ds
    real(kind = 8), dimension(neqn), intent(inout) :: y
    integer(kind = 4), intent(out) :: iflag

    ! local variables
    integer(kind = 4) :: i, ncount, n, idum, j
    real(kind = 8) :: x, htry

    real(kind = 8), dimension(neqn) :: yscal
    real(kind = 8), dimension(neqn) :: dydx
    real(kind = 8), dimension(neqn) :: yy

    real(kind = 8) ::hdid, hnext
    real(kind = 8) totdid, x1, x2, xx, hmin, tdum

    logical :: is_ok

    ! -------------------- executable --------------------
    ! store current values
    x=s	        !starting point
    htry=ds	!stepsize wanted to take

    ! interval to step through
    x1 = s        ! start of step
    x2 = s+ds	! end of step
    hmin=ds/maxstep	!minimum step size

    !save original values
    yy = y

    ! try rational function extrapolation first, we want to go from s to s+ds
    ncount=0	!# of times hdid < htry
    totdid=0.	!total stepsizes taken

    ! set initial error
    iflag=0

    is_ok = .false.

    if(err_tol.le.0.) then
       iflag=6
       write (6, *) ' Error tolerance eps = ', err_tol
       write (6, *)' <BS_ODE> Error tolerance <= 0.0, unable to continue'
    end if

    ! loop to calculaye the requested step
    do n=1,maxstep
       ! check if we overshot end pt, i.e. new x has to be within x1 and x2
       xx = x+htry
       if( (xx-x2)*(xx-x1) .gt. 0.0) htry=x2-x ! adjust htry to make sure the step stays withing x1 and x2

       ! get derivatives, dydx
       call derivs(tdum, yy, dydx, neqn)

       ! set up error scaling vector
       yscal=abs(yy) + abs(htry*dydx) + 1.e-30

       ! try to take a step
       call bsstep(yy,dydx,neqn,x,htry,yscal,hdid,hnext,derivs)  ! this also returns the new x-value
       ! call rk45_step(yy,dydx,neqn,x,htry,eps,yscal,hdid,hnext,derivs)

       if( (x-x2)*(x2-s).ge.0.) then
          is_ok = .true.
          exit	! on the limit or outside x1, x2 interval ------->quit loop
       endif

       ! x is where we have moved to in hdid stepsize
       if(hdid.lt.htry) then
          ! if it was not possible to makd a htry step with the required accucacy try a smaller step size fot bsstep
          ! mx number of bsstep tries 
          if(ncount.lt.max_bsstep_tries) then
             htry=hnext	!use the next recommended B-S stepsize
             cycle	!loop again
          else
             ! if bsstep failed to reach htry, then try a Runge-Kutta step
             ! starting at where bsstep left off
             htry=hnext/16.	!reduce step size for R-K
             do j=1,maxrk
                call derivs(tdum, yy, dydx, neqn)    !yy is changed, so need new dydx
                ! get new scaling variables              
                yscal = abs(yy) + abs(htry*dydx)+1.e-30	!so no 0's
                ! take a Runge_Kutta step
                call rk45_step(yy,dydx,neqn,x,htry,yscal,hdid,hnext,derivs)
                htry=hnext
             end do	!end runge kutta loop
             ncount=0
          end if	!if first bsstep loop or try runge kutta
       end if	!if didn't reach htry
       ! x and yy should be new values, hdid is how far we went
       ! now reset htry, loop back and try bsstep again
       htry=hnext
    end do	!end of maxsteps loop

    if (is_ok) then
       ! return new s, y
       s=x
       y = yy
    else
       ! if we get here, we didn't converge
       write (6, *)' bs_ode: could not converge on step size ', ds,' at this X location',x
       iflag = 5
    endif
    return
  end subroutine bs_ode


  subroutine bsstep(y,dydx,nv,x,htry,yscal,hdid,hnext,derivs)
    !     Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust step size.

    !     Input:
    !     y     -  dependent variable vector of length - nv
    !     dydx  -  derivatives at the staring value  of independent variable x
    !     nv    -  length of y and dydx
    !     htry  -  step size to be tried
    !     yscal -  scale to evaluating error     

    !     Output:
    !     y     - old and new values of y
    !     x     - old and new value of x
    !     hdid  - performed step size
    !     hnext - estimated next step size

    !     input function:
    !     derivs - user supplied function to calculate the derivatives      


    interface
       subroutine derivs(xx, yy, dyydxx, nn)
         integer(kind = 4), intent(in) :: nn
         real(kind = 8), intent(in) :: xx
         real(kind = 8), dimension(nn), intent(in) :: yy
         real(kind = 8), dimension(nn), intent(out) :: dyydxx
       end subroutine derivs
    end interface

    ! arguments
    real(kind = 8), dimension(nv), intent(inout) :: y
    real(kind = 8), dimension(nv), intent(in) ::dydx
    integer(kind = 4), intent(in) :: nv
    real(kind = 8), intent(inout) :: x
    real(kind = 8), intent(in) :: htry
    real(kind = 8), dimension(nv), intent(in) ::yscal
    real(kind = 8), intent(out) :: hdid
    real(kind = 8), intent(out)  :: hnext


    real(kind = 8), dimension(nv) ::yerr
    real(kind = 8), dimension(nv) ::ysav
    real(kind = 8), dimension(nv) ::dysav
    real(kind = 8), dimension(nv) ::yseq
    integer(kind = 4), dimension(imax) :: nseq = (/2,4,6,8,12,16,24,32,48,64,96/)
                                                 
    ! locals
    real(kind = 8) :: h, xsav, xest, scaled_err
    integer(kind = 4) :: i,j,k

    h=htry
    ! save initial values
    xsav = x
    ysav = y
    dysav = dydx

    do ! endless loop
       ! decrease step size
       do i=1,imax
          call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq,derivs)
          xest=(h/nseq(i))**2

          call rzextr(i,xest,yseq,y,yerr,nv,nuse)
          ! get the largest deviation
          scaled_err = maxval( abs(yerr/yscal) ) / err_tol

          if(scaled_err .lt. one) then
             x=x+h
             hdid=h
             if(i.eq.nuse)then
                hnext=h*shrink
             else if(i.eq.nuse-1)then
                hnext=h*grow
             else
                hnext=(h*nseq(nuse-1))/nseq(i)
             endif
             return ! normal case
          endif
       enddo
       h=0.25*h/2**((imax-nuse)/2)
       if(x+h.eq.x) pause 'step size underflow.'
    enddo

  end subroutine bsstep

  subroutine rzextr(iest,xest,yest,yz,dy,nv,n_use)

    !     Use diagonal rational function extrapolation (numerator and denominator polynomials with
    !     the same order) to evaluate nv functions at x = 0 by fitting a diagnonal rational function to
    !     a sequance of estimates with progressively smaller values x = xest and corresponding function
    !     vectors y = yest. This call is number iest in the sequence of calls. The extrapolation uses
    !     at most the last n_use estimates.
    !     Extrapolated function values are output as yz and the estimated error is output as dy


    integer(kind = 4), intent(in) :: iest
    real(kind = 8), intent(in) :: xest
    real(kind = 8), dimension(nv), intent(in) :: yest
    real(kind = 8), dimension(nv), intent(out) :: yz
    real(kind = 8), dimension(nv), intent(out) :: dy
    integer(kind = 4), intent(in) :: nv
    integer(kind = 4), intent(in) :: n_use

    ! locals 
    real(kind = 8), dimension(imax) :: x
    real(kind = 8), dimension(nmax,nuse) :: d
    real(kind = 8), dimension(nuse) :: fx
    real(kind = 8) :: b, b1, c, v, ddy, yy

    integer(kind = 4) :: i, j, k, m1

    x(iest)=xest
    if(iest.eq.1) then
       do j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
       enddo
    else
       m1=min(iest,n_use)
       do k=1,m1-1
          fx(k+1)=x(iest-k)/xest
       enddo
       do  j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do  k=2,m1
             b1=fx(k)*v
             b=b1-c
             if(b.ne.0.) then
                b=(c-v)/b
                ddy=c*b
                c=b1*b
             else
                ddy=v
             endif
             if(k.ne.m1) v=d(j,k)
             d(j,k)=ddy
             yy=yy+ddy
          enddo
          dy(j)=ddy
          yz(j)=yy
       enddo
    endif
    return
  end subroutine rzextr

  subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
    ! Modified midpoint step.
    ! Input:
    !    y: dependent variable vector of length nvar, evaluated at xs
    ! dydx: derivative of y at each value
    ! nvar: length of y, dydx
    !   xs: value of independent variable
    ! htot: total step to be made
    ! nstep: number of sub-steps to be used

    ! Output:
    ! yout: value of dependent varable vector at new location
    ! derivs: dyoutdx at the new location

    interface
       subroutine derivs(xx, yy, dyydxx, nn)
         integer(kind = 4), intent(in) :: nn
         real(kind = 8), intent(in) :: xx
         real(kind = 8), dimension(nn), intent(in) :: yy
         real(kind = 8), dimension(nn), intent(out) :: dyydxx
       end subroutine derivs
    end interface

    integer(kind = 4), intent(in) :: nvar
    integer(kind = 4), intent(in) :: nstep
    real(kind = 8), dimension(nvar), intent(in) :: y
    real(kind = 8), dimension(nvar), intent(in) :: dydx
    real(kind = 8), intent(in) :: xs
    real(kind = 8) :: htot
    real(kind = 8), dimension(nvar), intent(out) :: yout

    ! locals
    real(kind = 8), dimension(nvar) :: ym
    real(kind = 8), dimension(nvar) :: yn
    real(kind = 8), dimension(nvar) :: swap
    real(kind = 8) h, h2,  x

    integer(kind = 4) i, n

    h=htot/nstep
    ! operate on entire arrays
    ym = y
    yn = y + h*dydx
    x=xs+h
    call derivs(x, yn, yout, nvar)
    h2=2.*h
    do n=2,nstep
       swap = ym + h2*yout
       ym = yn
       yn = swap
       x=x+h
       call derivs(x,yn,yout, nvar)
    enddo
    yout = 0.5*(ym + yn + h*yout)
    return
  end subroutine mmid


  ! Fifth order Runge-Kutta step with monitoring of local truncation error
  ! to ensure accuracy and adjust stepsize.  

  ! Input:
  !	y	- dependent variable vector
  !	n	- size of vector
  !	dydx	- derivatives of Y
  !	x	- independent variable
  !	htry	- step size to be attempted (scaled down by 16 or 32 as nec.)
  !	eps	- required accuracy
  !	yscal	- vector of error scaling factors (|y(i)|+|h*dydx(i)|+1.e-30)
  !	derivs 	- user supplied subroutine to compute right-hand side derivs
  !
  ! output:
  !	y,x	- replaced with new values
  !	hdid	- stepsize actually accomplished
  !	hnext	- estimated next stepsize


  subroutine rk45_step(y,dydx,n,x,htry,yscal,hdid,hnext,derivs)

    interface
       subroutine derivs(xx, yy, dyydxx, nn)
         integer(kind = 4), intent(in) :: nn
         real(kind = 8), intent(in) :: xx
         real(kind = 8), dimension(nn), intent(in) :: yy
         real(kind = 8), dimension(nn), intent(out) :: dyydxx
       end subroutine derivs
    end interface

    ! arguments
    real(kind = 8), dimension(n), intent(inout) :: y
    real(kind = 8), dimension(n), intent(inout) :: dydx
    integer(kind = 4), intent(in) :: n
    real(kind = 8), intent(inout) :: x
    real(kind = 8), intent(in) :: htry
    real(kind = 8), dimension(n), intent(in) :: yscal
    real(kind = 8), intent(out) :: hdid
    real(kind = 8), intent(out) :: hnext

    ! local variables
    ! parameters
    integer(kind = 4), parameter :: nmax = 6  ! set for 6 dim. phase space
    real(kind = 8), parameter :: pgrow = -.20
    real(kind = 8), parameter :: pshrnk = -.25
    real(kind = 8), parameter :: fcor = 1./15.
    real(kind = 8), parameter :: one = 1.d0
    real(kind = 8), parameter :: safety = 0.9
    real(kind = 8), parameter :: errcon = 6.e-4 

    integer(kind = 4) :: istop
    integer(kind = 4) :: i
    real(kind = 8) :: errmax, h, hh, xdum
    real(kind = 8) :: xsav
    real(kind = 8), dimension(n) :: ytemp
    real(kind = 8), dimension(n) :: ysav
    real(kind = 8), dimension(n) :: dysav

    ! note:  errcon=(4/safety)**(1/pgrow)

    !-----------------------executable-----------------------

    ! save initial values
    xsav=x
    ysav=y
    dysav=dydx

    ! set stepsize to initial trial value
    h=htry

    do 
       ! take two half steps
       hh=0.5*h
       call rk45(ysav,dysav,n,xsav,hh,ytemp,derivs)
       
       x=xsav+hh ! advance x one half step
       call derivs(xdum,ytemp,dydx,n)
       call rk45(ytemp,dydx,n,x,hh,y,derivs)

       x=xsav+h ! advance x one half step
       if(x.eq.xsav) then
          write (6, *)' <rkqc> stepsize not significant'
       end if

       ! take the large step
       call rk45(ysav,dysav,n,xsav,h,ytemp,derivs)

       ! evaluate accuracy
       errmax=0.d0
       ytemp = y - ytemp	! error estimate
       errmax=max(errmax, maxval( abs(ytemp/yscal)) )
       errmax=errmax/err_tol	!scale relative to required tolerance
       if(errmax.gt.one) then	!truncation error too large, reduce stepsize
          h=safety*h*(errmax**pshrnk)
          cycle	!try again
       else	!step succeeded, compute next stepsize
          exit
       endif
    enddo
    hdid=h
    if(errmax.gt.errcon) then
       hnext=safety*h*(errmax**pgrow)
    else
       hnext=4.*h
    end if

    ! mop up fifth-order truncation error
    y = y + ytemp * fcor

    return
  end subroutine rk45_step

  !=======================================================
  ! Input:
  !	y	- dependent variable vector
  !	dydx	- derivatives of y
  !	n	- size of vector
  !	x	- independent variable
  !	h	- interval size
  !       derivs  - subroutine to calculate the derivatives (see interface)
  !
  ! Output:
  !       y_new   - dependent values at x + h 
  !
  !
  ! fourth-order Runge-Kutta method to advance the solution over an
  ! interval h. y_new can be the same as y

  subroutine rk45(y,dydx,n,x,h,y_new,derivs)

    interface
       subroutine derivs(xx, yy, dyydxx, nn)
         integer(kind = 4), intent(in) :: nn
         real(kind = 8), intent(in) :: xx
         real(kind = 8), dimension(nn), intent(in) :: yy
         real(kind = 8), dimension(nn), intent(out) :: dyydxx
       end subroutine derivs
    end interface

    integer(kind = 4) :: n
    real(kind = 8) :: x
    real(kind = 8) :: h
    real(kind = 8), dimension(n) :: y
    real(kind = 8), dimension(n) :: dydx
    real(kind = 8), dimension(n) :: y_new

    ! local variables
    integer(kind = 4) :: i

    real(kind = 8), parameter :: one_6th = 1.e0/6.
    real(kind = 8) :: hh, xh
    real(kind = 8), dimension(n) :: yt, dydx_t
    real(kind = 8), dimension(n) :: k1, k2, k3, k4

    !---------------------- executable ------------------
    ! check the array boundaries
    if (n .gt. nmax) then
       write (6,*) ' too many equations nmax = ', nmax, ' requested = ', n
       stop
    endif

    ! setup step sizes
    hh=h*0.5        ! half step
    
    xh=x+hh  ! new position 
    ! calculate first y increment from x = x + h/2, ! these are all arrays !
    k1 = h*dydx

    ! calculate derivatives for approx x + h/2 and y(x + h/2)
    call derivs(xh, y + 0.5*k1, dydx_t, n)
    k2 = h*dydx_t
    
    ! yt now evaluated at x + h/2 and y(x + h/2) estimated with dydx evaluated at x+h/2
    ! dyt now estimated at x + h/2 and y(x + h/2)

    ! third step: k3 
    call derivs(xh, y + 0.5*k2, dydx_t, n)
    k3 = h*dydx_t

    ! fourth step
    call derivs(x+h, y + k3, dydx_t, n)
    k4 = h*dydx_t

    ! accumulate increments with proper weights
    y_new = y + (k1 + 2.*k2 + 2.*k3 + k4)/6.d0

    return
  end subroutine rk45

end module bs_mod
