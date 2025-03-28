!> @brief Module containing the kind definitions for the real numbers
MODULE kinds

    IMPLICIT NONE

    !> Double precision kind
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)

    !> Single precision kind
    INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6)

END MODULE kinds
