program test_checklim

  use orbit_parameters_mod
  use flux_par_mod
  
  ! test that flux is read correctly and that check_lim works
  implicit none
  

  ! file name
  character(len = 132) :: efit_gfile_name = '029881.00252.dat'
  character(len = 132) :: efit_directory = '../example_data/'
  ! return values
  integer(kind = 4) :: i, j, k, get_flux, ier
  
  ! position vector
  real(kind = 8) :: rr, zz, phi
  real(kind = 8), dimension(1) :: xwant, ywant
  integer(kind = 4), dimension(nw*nh) :: is_inside

  ! field vector
  real(kind = 8) :: br, bz, bphi
  real(kind = 8) :: bmag

  ! get the flux data

  ier = get_flux(efit_directory, efit_gfile_name)
  print *, 'get_flux return value = ', ier
   
  
  ! test zlim over the entire flux grid
  call checklim (nw, nh, rgrid, zgrid, nbdry, rbdry, zbdry, is_inside)

  ! write results
  open(16, file = 'test_flux.data', status = 'unknown')
  write(16, *)'#! i[i,0]/ r[f,1]/ z[f, 2]/ inside[i,3]/'
  
  k = 0
  do i = 1, nw
     do j = 1, nh
        k = k + 1
        write(16, *) i, rgrid(i), zgrid(j), is_inside(k)
     enddo
  enddo
  close(16)
  
end program test_checklim

