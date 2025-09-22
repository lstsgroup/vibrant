MODULE pade

    USE gx_ac, ONLY: create_thiele_pade, evaluate_thiele_pade_at, params, free_params
    USE kinds, ONLY: dp
    IMPLICIT NONE

CONTAINS

    !> @brief interpolate n_parameter function values on a evenly
        !!        spaced grid of n_points using the thiele pade model
        !!        (only half of the given y values are used, zeros
        !!        in other half)
        !!
        !! @parameter[in]  n_parameter -- number of points y is tabulated
        !! @parameter[in]  y_ref       -- tabulated y values of a complex function
        !! @parameter[in]  n_points    -- number of interpolated points
        !! @parameter[out] y_out       -- tabulated interpolated values using pade
    SUBROUTINE interpolate(n_parameter, y_ref, n_points, y_out)
        INTEGER, INTENT(in)  :: n_parameter
        COMPLEX(kind=dp), DIMENSION(:), INTENT(in)  :: y_ref
        INTEGER, INTENT(in)  :: n_points
        COMPLEX(kind=dp), DIMENSION(:), INTENT(out) :: y_out

        ! internal variables
        TYPE(params) :: pade_params
        INTEGER      :: num_ref_points, i, n_halve_query_points
        INTEGER      :: last_important, last_important_out, n_par
        REAL(kind=dp) :: first, last, step, x_last_important
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x_ref_complx, x_ref_complx_tmp, y_ref_tmp
        COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x_out_complx
        ! convert x to complex type
        ! use only half of the array (other half filled with zeros)
        num_ref_points = SIZE(y_ref)
        ALLOCATE (x_ref_complx(num_ref_points))
        first = 1.0_dp
        last = 2.0_dp
        step = (last - first)/(num_ref_points - 1)
        DO i = 1, num_ref_points
            x_ref_complx(i) = COMPLEX(first + (i - 1.0_dp)*step, 0.0_dp)
        END DO
        last_important = num_ref_points
        DO i = num_ref_points, 1, -1
            IF (y_ref(i).EQ.COMPLEX(0.0_dp, 0.0_dp) .AND. i>INT(num_ref_points/2)) THEN
                !if (y_ref(i) .eq. complex(0.0_dp, 0.0_dp)) then
                last_important = i - 1
                x_last_important = REAL(x_ref_complx(i - 1), kind=dp)
            END IF
        END DO
        CALL select_points_evenly(x_ref_complx(1:last_important), last_important, 100d0, x_ref_complx_tmp, n_par)
        CALL select_points_evenly(y_ref(1:last_important), last_important, 100d0, y_ref_tmp, n_par)
        ! create model
        pade_params = create_thiele_pade(n_par, x_ref_complx_tmp, &
                                         y_ref_tmp, &
                                         do_greedy=.FALSE., PRECISION=64)

        ! create points where function is interpolated
        ALLOCATE (x_out_complx(n_points))
        step = (last - first)/(n_points - 1.0_dp)
        last_important_out = 1
        DO i = 1, n_points
            x_out_complx(i) = COMPLEX(first + (i - 1.0_dp)*step, 0.0_dp)
            IF (x_out_complx(i)%re.LE.x_last_important) THEN
                last_important_out = i
            END IF
        END DO

        ! evaluate model at given x for half the points (other half zero)
        y_out = evaluate_thiele_pade_at(pade_params, x_out_complx)
        y_out(last_important_out:) = COMPLEX(0.0_dp, 0.0_dp)

        ! deallocation
        CALL free_params(pade_params)
        DEALLOCATE (x_ref_complx)
        DEALLOCATE (x_out_complx)
    END SUBROUTINE interpolate
!> @brief Selects evenly spaced points from `points` according to percentage p
!> @param[in] points Input array of points (real numbers)
!> @param[in] N Total number of points in the input array
!> @param[in] p Percentage of points to select (0 < p <= 100)
!> @param[out] selected Output array of selected points (must be allocated outside)
!> @param[out] M Number of selected points (output)
    SUBROUTINE select_points_evenly(points, N, p, selected, M)
        INTEGER, INTENT(in) :: N
        COMPLEX(kind=8), DIMENSION(N), INTENT(in) :: points
        REAL(kind=8), INTENT(in) :: p
        COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE, INTENT(out) :: selected
        INTEGER, INTENT(out) :: M
        ! internal variables
        INTEGER :: i, step
        IF (p<=0.0d0 .OR. p>100.0d0) THEN
            PRINT *, "Error: percentage must be between 0 and 100."
            M = 0
            RETURN
        END IF
        ! Compute how many points to select
        M = INT((p/100.0d0)*DBLE(N) + 0.5d0)  ! rounded to nearest integer
        IF (M<1) M = 1
        IF (M>N) M = N
        ALLOCATE (selected(M))
        ! Select evenly spaced indices
        DO i = 1, M
            step = INT(DBLE(i - 1)*DBLE(N - 1)/DBLE(M - 1) + 0.5d0) + 1
            selected(i) = points(step)
        END DO
    END SUBROUTINE select_points_evenly
END MODULE pade
