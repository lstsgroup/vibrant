!> @brief proviede information on compilation / system setting
!!        used in a run
MODULE config_info

    USE omp_lib, ONLY: omp_get_max_threads

    IMPLICIT NONE

    PUBLIC :: output_config_info
    PRIVATE

    ! passing compiler settings from cmake to fortran
    INCLUDE "cmake_info.f90"

CONTAINS

    !> @brief report system/compilation configuration
    SUBROUTINE output_config_info()
        CALL output_vibrant_header()
        PRINT *, REPEAT('-', 30)
        CALL output_datetime()
        CALL output_used_threads()
        CALL output_compiler_settings()
        PRINT *, REPEAT('-', 30)
        PRINT *, ""
    END SUBROUTINE output_config_info

    !> @brief print the program name
    SUBROUTINE output_vibrant_header()
        PRINT *, ""
        PRINT *, " .--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--. "
        PRINT *, "/ .. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \"
        PRINT *, "\ \/\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ \/ /"
        PRINT *, " \/ /`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'\/ / "
        PRINT *, " / /\    _     _     __      _____     __ __       _____      __   __     _______    / /\ "
        PRINT *, "/ /\ \  /_/\ /\_\   /\_\   /\  __/\   /_/\__/\    /\___/\    /_/\ /\_\  /\_______)\ / /\ \"
        PRINT *, "\ \/ /  ) ) ) ( (   \/_/   ) )(_ ) )  ) ) ) ) )  / / _ \ \   ) ) \ ( (  \(___  __\/ \ \/ /"
        PRINT *, " \/ /  /_/ / \ \_\   /\_\ / / __/ /  /_/ /_/_/   \ \(_)/ /  /_/   \ \_\   / / /      \/ / "
        PRINT *, " / /\  \ \ \_/ / /  / / / \ \  _\ \  \ \ \ \ \   / / _ \ \  \ \ \   / /  ( ( (       / /\ "
        PRINT *, "/ /\ \  \ \   / /  ( (_(   ) )(__) )  )_) ) \ \ ( (_( )_) )  )_) \ (_(    \ \ \     / /\ \"
        PRINT *, "\ \/ /   \_\_/_/    \/_/   \/____\/   \_\/ \_\/  \/_/ \_\/   \_\/ \/_/    /_/_/     \ \/ /"
        PRINT *, " \/ /                                                                                \/ / "
        PRINT *, " / /\.--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--..--./ /\ "
        PRINT *, "/ /\ \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \.. \/\ \"
        PRINT *, "\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `'\ `' /"
        PRINT *, " `--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--'`--' "
        PRINT *, ""
    END SUBROUTINE output_vibrant_header

    !> @brief report how many OMP threads are available
    SUBROUTINE output_used_threads()
        INTEGER :: num_threads
        num_threads = omp_get_max_threads()
        PRINT *, "Number of OMP threads used: ", num_threads
    END SUBROUTINE output_used_threads

    !> @brief report current date and time
    SUBROUTINE output_datetime()
        INTEGER, DIMENSION(8) :: values
        CALL DATE_AND_TIME(values=values)
        PRINT '(A,I2.2,A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2)', &
            ' Date: ', values(3), '.', values(2), '.', values(1), &
            ', Time: ', values(5), ':', values(6), ':', values(7)
    END SUBROUTINE output_datetime

    !> @brief report compilation info
    SUBROUTINE output_compiler_settings()
        PRINT *, "Fortran compiler: ", compiler
        PRINT *, "Compiler flags: ", compiler_flags
        PRINT *, "FFTW library path: ", fft_lib_dir
        PRINT *, "GreenX library path: ", greenx_lib_dir
    END SUBROUTINE output_compiler_settings

END MODULE config_info
