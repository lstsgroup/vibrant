MODULE cell_types

    USE kinds, ONLY: dp
    USE constants, ONLY: pi, speed_light
    USE vib_types, ONLY: global_settings, systems

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: pbc, invert3x3, build_hmat, determinant3x3

CONTAINS

!*********************************************************************************************
!*********************************************************************************************

    SUBROUTINE build_hmat(sys, hmat)

        TYPE(systems), INTENT(IN) :: sys
        REAL(dp), INTENT(OUT) :: hmat(3, 3)
        REAL(dp) :: a, b, c, alpha, beta, gamma
        REAL(dp) :: ca, cb, cg, sa, sb, sg

        a = sys%cell%box_x
        b = sys%cell%box_y
        c = sys%cell%box_z

        alpha = sys%cell%angle_alpha*pi/180.0_dp
        beta = sys%cell%angle_beta*pi/180.0_dp
        gamma = sys%cell%angle_gamma*pi/180.0_dp

        ca = COS(alpha); cb = COS(beta); cg = COS(gamma)
        sa = SIN(alpha); sb = SIN(beta); sg = SIN(gamma)

        hmat(:, 1) = (/a, 0.0_dp, 0.0_dp/)
        hmat(:, 2) = (/b*cg, b*sg, 0.0_dp/)
        hmat(:, 3) = (/c*cb, c*(ca - cb*cg)/sg, &
                       c*SQRT(1.0_dp + 2.0_dp*ca*cb*cg - ca**2 - cb**2 - cg**2)/sg/)

    END SUBROUTINE build_hmat

!*********************************************************************************************
!*********************************************************************************************

    FUNCTION determinant3x3(a) RESULT(det)

        USE kinds, ONLY: dp
        REAL(dp), INTENT(IN) :: a(3, 3)
        REAL(dp) :: det

        det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) - &
              a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) + &
              a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))

    END FUNCTION determinant3x3

!*********************************************************************************************
!*********************************************************************************************

    SUBROUTINE invert3x3(a, ainv)

        USE kinds, ONLY: dp
        REAL(dp), INTENT(IN)  :: a(3, 3)
        REAL(dp), INTENT(OUT) :: ainv(3, 3)
        REAL(dp) :: det

        det = determinant3x3(a)

        IF (ABS(det)<1e-12_dp) STOP "Singular matrix in invert3x3"
        det = 1.0_dp/det

        ainv(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))*det
        ainv(1, 2) = (a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3))*det
        ainv(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))*det

        ainv(2, 1) = (a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3))*det
        ainv(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))*det
        ainv(2, 3) = (a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3))*det

        ainv(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))*det
        ainv(3, 2) = (a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2))*det
        ainv(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))*det

    END SUBROUTINE invert3x3

!*********************************************************************************************
!*********************************************************************************************

    SUBROUTINE pbc(coord2, coord1, sys, dr)

        USE kinds, ONLY: dp
        TYPE(systems), INTENT(INOUT) :: sys
        REAL(dp), DIMENSION(3), INTENT(IN)  :: coord2, coord1
        REAL(dp), DIMENSION(3), INTENT(OUT) :: dr
        REAL(dp) :: hmat(3, 3), h_inv(3, 3)
        REAL(dp) :: s(3), vec(3)

        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        sys%cell%vec = coord2 - coord1
        s = MATMUL(h_inv, sys%cell%vec)

        ! Minimum-image in fractional space ([-0.5,0.5))
        s = s - ANINT(s)

        dr = MATMUL(hmat, s)

    END SUBROUTINE pbc

END MODULE cell_types
