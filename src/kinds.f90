!> @brief Module containing the kind definitions for the real numbers
MODULE kinds

    IMPLICIT NONE

    !> Double precision kind
    INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(15)

    !> Single precision kind
    INTEGER, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND(6)

    !> default string length
    INTEGER, PARAMETER, PUBLIC :: str_len = 100

END MODULE kinds
