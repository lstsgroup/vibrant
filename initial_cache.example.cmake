# optional build configuration file for CMake
# uncommand non-default settings
# usage:
#       mkdir build && cd build
#       cmake -C initial_cache.example.cmake ..

# specifying the fortran compiler
#   set(CMAKE_Fortran_COMPILER "gfortran" CACHE FILEPATH "Fortran compiler")

# fortran compiler flags
#   set(CMAKE_Fortran_FLAGS "-O2 -fopenmp -ffree-line-length-none -cpp -I${CMAKE_BINARY_DIR}" CACHE STRING "Fortran compiler flags")

# disable regression tests
#   set(ENABLE_REGRESSION_TESTS OFF CACHE BOOL "Enable vibrant regression tests")

# Paths for external dependencies (optional)
#   set(BLAS_LIBRARY_DIR "/path/to/BLAS" CACHE PATH "Path to BLAS library directory")
#   set(LAPACK_LIBRARY_DIR "/path/to/LAPACK" CACHE PATH "Path to LAPACK library directory")
#   set(FFTW_INCLUDE_DIR "/path/to/fftw/install/include" CACHE PATH "Path to FFTW include directory")
#   set(FFTW_LIBRARY_DIR "/path/to/fftw/install/lib" CACHE PATH "Path to FFTW library directory")
#   set(GREENX_INCLUDE_DIR "/path/to/greenx/install/include/modules" CACHE PATH "Path to GreenX include directory")
#   set(GREENX_LIBRARY_DIR "/path/to/greenx/install/lib" CACHE PATH "Path to GreenX library directory")