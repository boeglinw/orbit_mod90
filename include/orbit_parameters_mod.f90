! module for orbit parameters and arrays
module orbit_parameters_mod

    implicit none
! old eparmdu129.h
      integer(kind = 4), parameter :: nfcoil=18
      integer(kind = 4), parameter :: nsilop=41
      integer(kind = 4), parameter :: nrogow=1
      integer(kind = 4), parameter :: ntime=1
      integer(kind = 4), parameter :: nacoil=1

      integer(kind = 4), parameter :: nesum=6
      integer(kind = 4), parameter :: magpri67=29
      integer(kind = 4), parameter :: magpri322=31
      integer(kind = 4), parameter :: magprit=6
      integer(kind = 4), parameter :: magpri=magpri67+magpri322
      
      integer(kind = 4), parameter :: necoil=122
      integer(kind = 4), parameter :: mpress=132
      integer(kind = 4), parameter :: nvesel=24
      integer(kind = 4), parameter :: ndata=41
      integer(kind = 4), parameter :: nwwcur=18
      integer(kind = 4), parameter :: nstark=16
      integer(kind = 4), parameter :: nffcur=18
      integer(kind = 4), parameter :: nppcur=18
      
      integer(kind = 4), parameter :: npcurn=nffcur+nppcur
      integer(kind = 4), parameter :: mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil
      integer(kind = 4), parameter :: npcur2=npcurn*2
      integer(kind = 4), parameter :: nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+mpress+nfcoil+nstark
      integer(kind = 4), parameter :: nwcurn=nwwcur+npcurn
      integer(kind = 4), parameter :: npcur3=npcurn*2
      integer(kind = 4), parameter :: nwcur2=nwcurn*2
     
      integer(kind = 4), parameter :: npoint=800
      integer(kind = 4), parameter :: nw=300
      integer(kind = 4), parameter :: nh=300
      
      integer(kind = 4), parameter :: nwnh=nw*nh
      integer(kind = 4), parameter :: nh2=2*nh
      integer(kind = 4), parameter :: nwrk=2*(nw+1)*nh
      integer(kind = 4), parameter :: ncurrt=nvesel+nesum+nfcoil
      
      integer(kind = 4), parameter :: mbdry=1500
      integer(kind = 4), parameter :: nbwork=nsilop
      integer(kind = 4), parameter :: kxiter=250
      integer(kind = 4), parameter :: mqwant=30
      integer(kind = 4), parameter :: nlimit=150
      integer(kind = 4), parameter :: nlimbd=6
      
      integer(kind = 4), parameter :: msbdry=mbdry+nsilop+nfcoil+1
      integer(kind = 4), parameter :: msbdr2=2*msbdry
      integer(kind = 4), parameter :: nrsma2=2*nrsmat
      integer(kind = 4), parameter :: nwwf=2*nw
      integer(kind = 4), parameter :: nwf=nwwf
      
      integer(kind = 4), parameter :: nxtram=10
      integer(kind = 4), parameter :: nxtrap=npoint
      integer(kind = 4), parameter :: nxtlim=9
      integer(kind = 4), parameter :: nco2v=3
      integer(kind = 4), parameter :: nco2r=1
      integer(kind = 4), parameter :: nangle=64
      integer(kind = 4), parameter :: nfbcoil=12
      integer(kind = 4), parameter :: kubicx = 4
      integer(kind = 4), parameter :: kubicy = 4
      integer(kind = 4), parameter :: lubicx = nw - kubicx + 1
      integer(kind = 4), parameter :: lubicy = nh - kubicy + 1
      integer(kind = 4), parameter :: kujunk = kubicx*kubicy*lubicx*lubicy
      
      integer(kind = 4), parameter :: modef=4
      integer(kind = 4), parameter :: modep=4
      integer(kind = 4), parameter :: modew=4 
      integer(kind = 4), parameter :: kubics=4 
      integer(kind = 4), parameter :: nrsp=200

      integer(kind = 4), parameter :: nxknot = lubicx+2*kubicx-1
      integer(kind = 4), parameter :: nyknot = lubicy+2*kubicy-1

      ! max. number of detectors
      integer(kind = 4), parameter :: n_par = 20

      ! some control parameters
      integer(kind = 4), parameter :: n111 = 1
      integer(kind = 4), parameter :: n333 = 3

      ! io unit numbers
      integer(kind = 4), parameter:: neqdsk = 38
      
      ! constants
      real(kind = 8), parameter :: pi  = 3.1415926535897932d0
      real(kind = 8), parameter :: twopi=2.0*pi
      real(kind = 8), parameter :: dtr = pi/180.

      ! conversion factor for magnetic moment
      real(kind = 8), parameter :: tmu = 2.0e-07
      real(kind = 8) XXXXXX  ! DUMMY valriable
end module orbit_parameters_mod
