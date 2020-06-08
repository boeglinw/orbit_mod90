
function get_flux(idir, ifname)
	  
  !**********************************************************************
  !                                                                  
  !     SUBPROGRAM DESCRIPTION:                                      
  !         - composes the full path for the eqdsk file and reads it
  !         - also add the g for the standard file name convention
  !     RECORD OF MODIFICATION:                                      
  !          
  !          06/2020 new version for f90 W. Boeglin                              
  !                                                                  
  !**********************************************************************


  use orbit_parameters_mod

  use flux_par_mod
  
  implicit none
  ! input variables
  character(len = *), intent(in):: idir
  character(len = *), intent(in):: ifname

  ! local variables
  character(len = 132):: std_fname
  integer(kind=4):: ire
  
  ! control variables for rdpar
  integer(kind=4):: imfit

  integer(kind = 4) :: get_flux

  ! for testing to make sure name are correct
	       
  print *, '---> get_flux  ifname: ', ifname
  print *, '---> get_flux   ifdir: ', idir

  ! setup standard file name remove leading and trailing spaces
  std_fname= TRIM(ADJUSTL(idir))//'/g'//ifname
  
  print *, '---> get_flux EQDSK file : '//std_fname
  
  call read_eqdsk(std_fname,imfit,ire)

  get_flux = ire
  
  return

end function get_flux

