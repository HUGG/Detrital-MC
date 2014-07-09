module mod_constants

implicit none

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  real(kind=dp),parameter :: pi = 3.141592653589793238462643383279502884197_dp

end module mod_constants
