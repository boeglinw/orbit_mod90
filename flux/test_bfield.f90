program test_bfield
  use em_fields_mod

  ! test bfield
  implicit none
  
  integer(kind = 4) :: get_flux

  ! file name
  character(len = 132) :: efit_gfile_name = '029881.00252.dat'
  character(len = 132) :: efit_directory = '../example_data/'
  ! return values
  integer(kind = 4) :: imfit, ier, i
  
  ! position vector
  real(kind = 8) :: rr, zz, phi

  ! field vector
  real(kind = 8), dimension(5):: b_vect

   ! get the flux data

  ier = get_flux(efit_directory, efit_gfile_name)
  print *, 'get_flux return value = ', ier

  
  do i = 1, 5
     print *, 'enter position nr. # ', i, ' to get the field :'
     read (*,*) rr, zz
     print *, ' r = ', rr, ', z = ', zz
     b_vect = bfield(rr, zz)
     print *, 'field vector : ', b_vect(1), b_vect(2),  b_vect(3)
     print *, 'magnitude : ', b_vect(4)
     print *, 'relative flux : ', b_vect(5)
  enddo
     
end program test_bfield
