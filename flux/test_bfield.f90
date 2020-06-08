program test_bfield

  ! test bfield
  implicit none

  interface
     function  bfield (r, z)  result(b)
       real(kind = 8), intent(in) :: r, z
       real(kind = 8), dimension(4) :: b
     end function bfield
  end interface
  
  integer(kind = 4) :: get_flux

  ! file name
  character(len = 132) :: efit_gfile_name = '029881.00252.dat'
  character(len = 132) :: efit_directory = '../../orbit_public/MAST_efit/'
  ! return values
  integer(kind = 4) :: imfit, ier, i
  
  ! position vector
  real(kind = 8) :: rr, zz, phi

  ! field vector
  real(kind = 8), dimension(4):: b_vect

   ! get the flux data

  ier = get_flux(efit_directory, efit_gfile_name)
  print *, 'get_flux return value = ', ier

  
  do i = 1, 5
     print *, 'enter position nr. # ', i, ' to get the field :'
     read (*,*) rr, zz
     print *, ' r = ', rr, ', z = ', zz
     b_vect = bfield(rr, zz)
     print *, 'field vector : ', b_vect(1), b_vect(2),  b_vect(3)
     print *, ' magnitude : ', b_vect(4)
  enddo
     
end program test_bfield
