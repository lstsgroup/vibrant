!
!   Copyright 2025 Ekin E. Winogradow, Johannes Scheffler, Moritz Leucke, Dorothea Golze
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!

!> @brief Module containing all procedures involving periodic boundary conditions
MODULE cell_types

    USE kinds, ONLY: dp
    USE constants, ONLY: pi, speed_light
    USE vib_types, ONLY: global_settings, systems

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: pbc, pbc_vec, invert3x3, build_hmat, determinant3x3

CONTAINS

!*********************************************************************************************
!*********************************************************************************************

    !> @brief Constructs the 3×3 cell matrix (h-matrix) from lattice parameters.
    !>
    !> This subroutine builds the transformation matrix `hmat` that defines the
    !> simulation cell vectors in Cartesian coordinates, based on the lattice
    !> constants (a, b, c) and angles (alpha, beta, gamma) provided in `sys%cell`.
    !>
    !> @param[in]  sys  -- System data structure (provides cell parameters)
    !> @param[out] hmat -- 3×3 matrix of cell vectors in Cartesian coordinates
    !>
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

    !> @brief Computes the determinant of a 3×3 matrix.
    !>
    !> @param[in]  a   -- 3×3 real matrix
    !> @return     det -- Determinant of the input matrix
    !>
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

    !> @brief Computes the inverse of a 3×3 matrix.
    !>
    !> @param[in]  a    -- 3×3 real input matrix
    !> @param[out] ainv -- 3×3 real inverse of the input matrix
    !>
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

    !> @brief Computes the minimum-image displacement vector between two points under
    !>        periodic boundary conditions (PBC).
    !>
    !>   Applies the minimum image convention by mapping the Cartesian displacement
    !>   between `coord2` and `coord1` into fractional coordinates, wrapping it into
    !>   the primary simulation cell, and converting back to Cartesian space.
    !>
    !> @param[inout] sys     -- System structure (provides the cell parameters)
    !> @param[in]    coord2  -- Cartesian coordinates of atom 2
    !> @param[in]    coord1  -- Cartesian coordinates of atom 1
    !> @param[out]   dr      -- Minimum-image displacement vector (coord2 - coord1) under PBC
    !> 
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
!*********************************************************************************************
!*********************************************************************************************

    !> @brief Computes distances between a reference coordinate and multiple points
    !>        under periodic boundary conditions (PBC).
    !>
    !> Applies the minimum image convention for each point in `coord_vec` relative to
    !> `coord1`, ensuring that distances are measured within the simulation
    !> cell. Returns the scalar distances for all input points.
    !>
    !> @param[inout] sys        -- System structure (provides the cell parameters)
    !> @param[in]    n_points   -- Number of points to evaluate
    !> @param[in]    coord_vec  -- Cartesian coordinates of target points (3 × n_points)
    !> @param[in]    coord1     -- Reference Cartesian coordinate
    !> @param[out]   distances  -- Minimum-image scalar distances from `coord1` to each point
    !>
    SUBROUTINE pbc_vec(n_points, coord_vec, coord1, sys, distances)

        TYPE(systems), INTENT(INOUT) :: sys
        INTEGER, INTENT(in) :: n_points
        REAL(dp), DIMENSION(3, n_points), INTENT(IN)  :: coord_vec
        REAL(dp), DIMENSION(3), INTENT(IN)  :: coord1
        REAL(dp), DIMENSION(n_points), INTENT(OUT) :: distances

        INTEGER :: i_point
        REAL(dp) :: hmat(3, 3), h_inv(3, 3)
        REAL(dp) :: s(3), vec(3)
        REAL(dp), DIMENSION(3, n_points) :: dr_list
        
        CALL build_hmat(sys, hmat)
        CALL invert3x3(hmat, h_inv)

        DO i_point = 1, n_points

            vec = coord_vec(:, i_point) - coord1

            s(1) = h_inv(1,1)*vec(1) + h_inv(1,2)*vec(2) + h_inv(1,3)*vec(3)
            s(2) = h_inv(2,1)*vec(1) + h_inv(2,2)*vec(2) + h_inv(2,3)*vec(3)
            s(3) = h_inv(3,1)*vec(1) + h_inv(3,2)*vec(2) + h_inv(3,3)*vec(3)

            s = s - ANINT(s)

            dr_list(1, i_point) = hmat(1,1)*s(1) + hmat(1,2)*s(2) + hmat(1,3)*s(3)
            dr_list(2, i_point) = hmat(2,1)*s(1) + hmat(2,2)*s(2) + hmat(2,3)*s(3)
            dr_list(3, i_point) = hmat(3,1)*s(1) + hmat(3,2)*s(2) + hmat(3,3)*s(3)

            distances(i_point) = SQRT(dr_list(1, i_point)**2 + dr_list(2, i_point)**2&
                                      + dr_list(3, i_point)**2)
        END DO

    END SUBROUTINE pbc_vec

END MODULE cell_types
