!> @brief Module containing the kind definitions for the real numbers
module kinds

  implicit none

  !> Double precision kind
  integer, parameter :: dp = selected_real_kind(15)

  !> Single precision kind
  integer, parameter :: sp = selected_real_kind(6)

end module kinds