module control_mod
  ! these are control parameters for controlling the calculations in the various trackers
  logical :: time_reversed = .true.   ! do a time reversed calculation (inverses the B-field direction)
  logical :: check_limiter = .true.   ! check if the calculated position is witing a limiter
  logical :: print_hit = .false. ! if true print a message when a limiter is hit
end module control_mod
