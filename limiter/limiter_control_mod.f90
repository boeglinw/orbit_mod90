module limiter_control_mod
  ! this is an interface module to control limiter setting
  use limiter_mod

contains
  ! initialization and  wrapper functions
  subroutine set_print_hit
    print_hit = .True.
  end subroutine set_print_hit

  subroutine clear_print_hit
    print_hit = .False.
  end subroutine clear_print_hit
  
  ! limiter initialization routines
  subroutine set_limiter_file_name(file_name)
    implicit none
    character(len = *), intent(in) :: file_name
     limiter_filename = file_name
  end subroutine set_limiter_file_name

  subroutine set_limiter_directory(directory)
    implicit none
    character(len = *), intent(in) :: directory
    limiter_directory = directory
  end subroutine set_limiter_directory
  
  subroutine limiter_init
    ! wrapping function for init_limiter
    call init_limiter
    return
  end subroutine limiter_init

end module limiter_control_mod
