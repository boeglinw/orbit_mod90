subroutine read_eqdsk(fname,imfit,ier)
  
  !**********************************************************************
  !**                                                                  **
  !**                                                                  **
  !**     SUBPROGRAM DESCRIPTION:                                      **
  !**          read_eqdsk reads  out the GAQ type eqdsk.               **
  !**                                                                  **
  !**     CALLING ARGUMENTS:                                           **
  !**                                                                  **
  !**     REFERENCES:                                                  **
  !**          (1)                                                     **
  !**          (2)                                                     **
  !**                                                                  **
  !**     RECORD OF MODIFICATION:                                      **
  !**          29/06/83..........first created                         **
  !**          02/02/99..........modified by DSD for PPPL orbit code   **
  !**          06/21/2018.......WB made f90 version removed            **
  !**                           code for ioption != 1                  **
  !**                                                                  **
  !**********************************************************************
  
  use orbit_parameters_mod
  
  use flux_par_mod
  
  implicit none
  
  
  integer(kind=4), parameter :: npitch=8*nlimit

  ! arguments
  character(len=*), intent(in):: fname
  integer(kind=4), intent(out):: imfit
  integer(kind=4), intent(out):: ier

  ! local variables

  character(len=10), dimension(6):: case
  
  real(kind = 8), dimension(nw):: workk
  real(kind = 8), dimension(nw):: xsi
  real(kind = 8), dimension(6):: pds(6)
  
  real(kind = 8), dimension(nw,nh):: psirz
  
  real(kind = 8):: zmid
  real(kind = 8):: rdim, zdim
  real(kind = 8):: ssibry, ssimag
  
  real(kind = 8):: darea, drgrid, dxsi, dzgrid
   
  integer(kind = 4):: negcur
  integer(kind = 4):: ioerr


  ! dummy variables (to satify file format)
  real(kind = 8) :: xdum, ydum, sdum
  
  ! loop counters
  integer(kind = 4):: i, j, k, kk 

  !-----------------------------Executable code begins here-----------------
  
  print*,' Entered routine read_eqdsk which is called from flux_init'
  
  ier = 0
  
  
  ! called from fluxinit then fnamein= fname2= g + fname= g+ well.= gwell.
  ! this fname is stored from scanned input file in subroutine rdpar
  print *, '----------------------------------------------------------------------'
  print *, '--> EQ File unit, name : ', neqdsk, fname
  open(unit=neqdsk,file=fname,status='old', iostat = ioerr)
  if (ioerr .ne. 0) then
     print *, 'problem opening ', fname, ' status = ', ioerr
     ier = 1
     return
  endif
  
  ! start reading the file
  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(6a8,3i4)') (case(i),i=1,6),imfit,mw,mh

  ! print information
  ! write (8, '('' Case='', 6a10)') (case(i), i=1,6)

  print*,'case(1) = ',case(1)
  print*,'case(2) = ',case(2)
  print*,'case(3) = ',case(3)
  print*,'case(4) = ',case(4)
  print*,'case(5) = ',case(5)
  print*,'case(6) = ',case(6)
  print*,'imfit = ',imfit
  print*,'mw = ',mw
  print*,'mh = ',mh

  ! case(i): identifier strings
  ! imfit : variable, controlling the re-scaling of ffprime amd pprime
  ! mw: number of hor. R-grid points
  ! mh: number of vertical. Z-grid points
  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(5e16.9)') rdim,zdim,rzero,rgrid(1),zmid

  print*,'rdim = ',rdim
  print*,'zdim = ',zdim
  print*,'rzero = ',rzero 
  print*,'rgrid(1) = ',rgrid(1)
  print*,'zmid = ',zmid

  ! xdim : hor. (R) dimension of computational box (m)
  ! zdim : ver. (Z) dimension of computational box (m)
  ! rzero: R in meter of vacuum toroidal magnetic field BCENTR
  ! rgrid(1) : Minimum R in meter of rectangular computational box
  ! zmid: Z of center of computational box in meter
  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(5e16.9)') rmaxis,zmaxis,ssimag,ssibry,bzero
  print*,'rmaxis = ',rmaxis 
  print*,'zmaxis = ',zmaxis 
  print*,'ssimag = ',ssimag 
  print*,'ssibry = ',ssibry 

  ! rmaxis: R of magnetic axis in meter
  ! zmaxis : Z of magnetic axis in meter
  ! ssimag : poloidal flux at magnetic axis in Weber/rad
  ! ssibry: poloidal flux at the plasma boundary in Weber /rad
  ! bzero: Vacuum toroidal magnetic field in Tesla at rzero

  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(5e16.9)') current,xdum,xdum,rmaxis,xdum
  print*,'current= ',current

  ! piii: Plasma current in Ampere

  ! all the other variables are basically read before
  ! xdum: dummy variables

  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(5e16.9)') zmaxis,xdum,sdum,xdum,xdum
  print *, 'zmaxis = ', zmaxis
  
  ! zmaxis: Z of magnetic axis in meter
  ! rest are dummy variables for parameters already read
  
  print *, '----------------------------------------------------------------------'
  print *, 'reading fpol '
  READ(neqdsk, '(5e16.9)') (fpol(i),i=1,mw)
  print*,'fpol(1)= ',fpol(1)

  ! fpol: Poloidal current function F in m-T, F = RB_T on flux grid
  
  print *, '----------------------------------------------------------------------'
  print *, 'reading pres '
  READ(neqdsk, '(5e16.9)') (pres(i),i=1,mw) 
  print*,'PRES(1)= ',pres(1) 

  ! pres: Plasma pressure in nt / m^2 on uniform flux grid
  
  print *, '----------------------------------------------------------------------'
  print*,'reading FFRIME into workk array'  
  READ(neqdsk, '(5e16.9)') (workk(i),i=1,mw) 
  print*,'workk(1)= ',workk(1)

  ! workk: FF’(ψ) in (mT)2 / (Weber /rad) on uniform flux grid
  !
  ! rescale ffprime according top the switch imfit
  ! note ffprime changes sign in any case
  !
  do i=1,mw
     if (imfit.ge.0) ffprim(i)=-workk(i)/(twopi*tmu)
     if (imfit.lt.0) ffprim(i)=-workk(i)
  enddo
  print*,'Finished reading FFRIME'  

  
  print *, '----------------------------------------------------------------------'
  print*,'reading PPRIME into workk array'
  READ(neqdsk, '(5e16.9)') (workk(i),i=1,mw)
  
  ! workk: P’(ψ) in (nt /m2) / (Weber /rad) on uniform flux grid
  ! change sign of pprime`
  do i=1,mw
     pprime(i)=-workk(i)
  enddo
  print*, 'pprime(1) = ', pprime(1)
  
  print*,'Finished reading PPRIME'  
      
  print *, '----------------------------------------------------------------------'
  print*, 'reading psirz'
  READ(neqdsk, '(5e16.9)') ((psirz(i,j),i=1,mw),j=1,mh)
  print*, 'finished reading psirz'


  print*, '----------------------------------------------------------------------'
  print*, 'reading QPSI' 
  READ(neqdsk, '(5e16.9)') (qpsi(i),i=1,mw)
  print*,'qpsi(1)= ', qpsi(1)
  
  ! qpsi:  q values on uniform flux grid from axis to boundary
  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(2i5)') nbdry,limitr
  print *, 'reading plasma boundary'
  print *, 'nbdry, limitr = ', nbdry,limitr
  
  ! nbdry: Number of boundary points
  ! limitr: Number of limiter points
  
  print *, '----------------------------------------------------------------------'
  READ(neqdsk, '(5e16.9)') (rbdry(i),zbdry(i),i=1,nbdry)

  ! rbdry: R of boundary points in meter
  ! zbdry: Z of boundary points in meter
  
  print *, '----------------------------------------------------------------------'
  print *, 'reading limiter data'
  READ(neqdsk, '(5e16.9)') (rlim(i),zlim(i),i=1,limitr)

  ! xlim: R of surrounding limiter contour in meter
  ! ylim: Z of surrounding limiter contour in meter

  print *, '----------------------------------------------------------------------'
  print *, ' ---> All read, setup data '
  print *, '----------------------------------------------------------------------'

  ! store local variables in final variables
  simag=ssimag
  psibry=ssibry
  sifm=simag
  sifb=psibry
  
  ! psirz: Poloidal flux in Weber / rad on the rectangular grid points
  ! store in a 1d array psi

  kk = 0
  do i=1,mw
     do j=1,mh
        kk=(i-1)*mh+j
        psi(kk)=psirz(i,j)
     enddo
  enddo
  print*, 'psi(1) = ', psi(1)

  ! setup the grid values
  drgrid=rdim/float(mw-1)
  dzgrid=zdim/float(mh-1)
  do i=1,mw
     rgrid(i)=rgrid(1)+(i-1)*drgrid
  enddo
  do i=1,mh
     zgrid(i)=zmid-zdim/2.+(i-1)*dzgrid
  enddo
  darea=drgrid*dzgrid
  print*,'rgrid(mw)= ',rgrid(mw)
  print*,'zgrid(1)= ',zgrid(1)
  print*,'zgrid(mh)= ',zgrid(mh)
  
  ! find limits
  
  rlmin=rlim(1)
  rlmax=rlmin
  zlmin=zlim(1)
  zlmax=zlmin
  do i=2,limitr
     rlmin=min(rlmin,rlim(i))
     rlmax=max(rlmax,rlim(i))
     zlmin=min(zlmin,zlim(i))
     zlmax=max(zlmax,zlim(i))
  enddo
  
  close(unit=neqdsk)

  ! setup relative variable Xsi from 0. to 1. for the radial grid 
  dxsi=1./float(mw-1)
  do i=1,mw
     xsi(i)=(i-1)*dxsi
     xxxsi(i)=xsi(i)
  enddo
  close(unit=neqdsk)
  

  !------------------------------------------------------------------------------
  ! end of reading EQDSK file
  !------------------------------------------------------------------------------
  
  !----------------------------------------------------------------------
  !    setup 2-d zpline for psi 
  !----------------------------------------------------------------------
  call sets2d(psi,c,rgrid,mw,bkx,lkx,zgrid,mh,bky,lky,wk,ier)

  !----------------------------------------------------------------------
  !   setup zpline for poliodal current function F
  !----------------------------------------------------------------------
  
  mwfpol=mw
  call zpline(mwfpol,xxxsi,fpol,bfpol,cfpol,dfpol)

  print *, '----------------------------------------------------------------------'
  print *, ' All DONE !'
  print *, '----------------------------------------------------------------------'

  
  return
  
end subroutine read_eqdsk

    
      
