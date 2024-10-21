module pade

    use gx_ac, only: create_thiele_pade, evaluate_thiele_pade_at, params, free_params

    implicit none

    contains 

        !> @brief interpolate n_parameter function values on a evenly 
        !!        spaced grid of n_points using the thiele pade model  
        !!        (only halve of the given y values are used, zeros 
        !!        in other halve)
        !!
        !! @parameter[in]  n_parameter -- number of points y is tabulated
        !! @parameter[in]  y_ref       -- tabulated y values of a complex function
        !! @parameter[in]  n_points    -- number of interpolated points 
        !! @parameter[out] y_out       -- tabulated interpolated values using pade
        subroutine interpolate(n_parameter,y_ref, n_points, y_out)
            integer,                       intent(in)  :: n_parameter
            complex(kind=8), dimension(:), intent(in)  :: y_ref
            integer,                       intent(in)  :: n_points
            complex(kind=8), dimension(:), intent(out) :: y_out

            ! internal variables
            type(params) :: pade_params
            integer      :: num_ref_points, i, n_halve_query_points
            integer      :: last_important, last_important_out
            real(kind=8) :: first, last, step, x_last_important
            complex(kind=8), dimension(:), allocatable :: x_ref_complx
            complex(kind=8), dimension(:), allocatable :: x_out_complx

            ! convert x to complex type 
            ! use only halve of the array (other halve filled with zeros)
            num_ref_points = size(y_ref)
            allocate(x_ref_complx(num_ref_points))
            first = 1.0d0
            last = 2.0d0
            step = (last - first) / (num_ref_points - 1)
            do i = 1, num_ref_points
                x_ref_complx(i) = complex(first + (i-1.0d0)*step, 0.0d0)
            end do
            last_important = num_ref_points
            do i = num_ref_points, 1, -1
                if (y_ref(i) .eq. complex(0.0d0, 0.0d0) .AND. i>int(num_ref_points/2)) then
                !if (y_ref(i) .eq. complex(0.0d0, 0.0d0)) then
                    last_important = i -1
                    x_last_important = real(x_ref_complx(i-1), kind=8)
                end if 
            end do

            ! create model
            pade_params = create_thiele_pade(last_important, x_ref_complx(1:last_important), &
                                             y_ref(1:last_important),      &
                                             do_greedy = .false., precision=64)

            ! create points where function is interpolated
            allocate(x_out_complx(n_points))
            step = (last - first) / (n_points - 1.0d0)
            last_important_out = 1
            do i = 1, n_points
                x_out_complx(i) = complex(first + (i-1.0d0)*step, 0.0d0)
                if (x_out_complx(i)%re .le. x_last_important) THEN
                    last_important_out = i
                end if 
            end do 

            ! evaluate model at given x for halve the points (other halve zero)
            y_out = evaluate_thiele_pade_at(pade_params, x_out_complx)
            y_out(last_important_out:) = complex(0.0d0, 0.0d0)

            ! deallocation 
            call free_params(pade_params)
            deallocate(x_ref_complx)
            deallocate(x_out_complx)

        end subroutine interpolate

end module pade
