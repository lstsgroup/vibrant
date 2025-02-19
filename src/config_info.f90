!> @brief proviede information on compilation / system setting 
!!        used in a run
module config_info

    use omp_lib, only: omp_get_max_threads

    implicit none 

    public :: output_config_info
    private 

    ! passing compiler settings from cmake to fortran
    include "cmake_info.f90"

    contains 

        !> @brief report system/compilation configuration
        subroutine output_config_info()
            call output_vibrant_header()
            print *, repeat('-', 30)
            call output_datetime()
            call output_used_threads()
            call output_compiler_settings()
            print *, repeat('-', 30)
            print *, ""
        end subroutine output_config_info

        !> @brief print the program name
        subroutine output_vibrant_header()
            print *, ""
            print *, " .--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--. "
            print *, "/ .. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \"
            print *, "\ \/\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ \/ /"
            print *, " \/ /`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'\/ / "
            print *, " / /\    _     _     __      _____     __ __       _____      __   __     _______    / /\ "
            print *, "/ /\ \  /_/\ /\_\   /\_\   /\  __/\   /_/\__/\    /\___/\    /_/\ /\_\  /\_______)\ / /\ \"
            print *, "\ \/ /  ) ) ) ( (   \/_/   ) )(_ ) )  ) ) ) ) )  / / _ \ \   ) ) \ ( (  \(___  __\/ \ \/ /"
            print *, " \/ /  /_/ / \ \_\   /\_\ / / __/ /  /_/ /_/_/   \ \(_)/ /  /_/   \ \_\   / / /      \/ / "
            print *, " / /\  \ \ \_/ / /  / / / \ \  _\ \  \ \ \ \ \   / / _ \ \  \ \ \   / /  ( ( (       / /\ "
            print *, "/ /\ \  \ \   / /  ( (_(   ) )(__) )  )_) ) \ \ ( (_( )_) )  )_) \ (_(    \ \ \     / /\ \"
            print *, "\ \/ /   \_\_/_/    \/_/   \/____\/   \_\/ \_\/  \/_/ \_\/   \_\/ \/_/    /_/_/     \ \/ /"
            print *, " \/ /                                                                                \/ / "
            print *, " / /\.--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--./ /\ "
            print *, "/ /\ \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \/\ \"
            print *, "\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `' /"
            print *, " `--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--' "
            print *, ""
        end subroutine output_vibrant_header

        !> @brief report how many OMP threads are available  
        subroutine output_used_threads()
            integer :: num_threads
            num_threads = omp_get_max_threads()
            print *, "Number of OMP threads used: ", num_threads             
        end subroutine output_used_threads

        !> @brief report current date and time
        subroutine output_datetime()
            integer, dimension(8) :: values
            call date_and_time(values=values)
            print '(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)', &
                ' Date: ', values(3), '.', values(2), '.', values(1), &
                ', Time: ', values(5), ':', values(6), ':', values(7) 
        end subroutine output_datetime
        
        !> @brief report compilation info
        subroutine output_compiler_settings()
            print *, "Fortran compiler: ", compiler 
            print *, "Compiler flags: ", compiler_flags
            print *, "FFTW library path: ", fft_lib_dir
            print *, "GreenX library path: ", greenx_lib_dir
        end subroutine output_compiler_settings

end module config_info