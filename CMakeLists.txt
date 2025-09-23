cmake_minimum_required(VERSION 3.18)

# set default fortran compiler to gfortran
if(NOT CMAKE_Fortran_COMPILER)
    message(STATUS "no fortran compiler specified")
    message(STATUS "using the default compiler: gfortran")
    set(CMAKE_Fortran_COMPILER gfortran)
else()
    message(STATUS "using specified fortran compiler: ${CMAKE_Fortran_COMPILER}")
endif()

# Project name and Fortran standard
project(vibrant Fortran)

# cmake module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Website 
option(ENABLE_DOCS "enable website generation" OFF)
message(STATUS "Enable website generation: ${ENABLE_DOCS}")
if (${ENABLE_DOCS})
    include(cmake/website.cmake)
    return()
endif()

# external dependencies
include(ExternalProject)
include(cmake/fftw_dependency.cmake)
include(cmake/greenx_dependency.cmake)
include(cmake/lapack_dependency.cmake)
set(LIBS fftw3 gx_ac ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

## Compiler flags
if(NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "-O -g -fopenmp -ffree-line-length-none -fcheck=bounds,do,recursion,pointer -cpp -I${CMAKE_BINARY_DIR}")
    message(STATUS "CMAKE_Fortran_FLAGS not defined. Using default: ${CMAKE_Fortran_FLAGS}")
else()
    message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif()

# Source files
set(SOURCES
    src/kinds.f90
    src/timing.f90
    src/constants.f90
    src/read_input.f90
    src/config_info.f90
    src/vib_types.f90
    src/setup.f90
    src/pade_interpolation.f90
    src/read_traj.f90
    src/dipole_calc.f90
    src/vel_cor.f90
    src/fin_diff.f90
    src/spectra.f90
    src/main.f90
)

# Executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} ${LIBS})

# Ensure FFTW and GreenX are built before linking
add_dependencies(${PROJECT_NAME} fftw greenx)

# regression tests
option(ENABLE_REGRESSION_TESTS "Enable vibrant regression tests" ON)
message(STATUS "Enable regression tests: ${ENABLE_REGRESSION_TESTS}")
if (${ENABLE_REGRESSION_TESTS})
    include(cmake/regression_tests.cmake)
endif()

# save cmake variables in f90 file to make them accesseble from fortran
configure_file(src/cmake_info.f90.in ${CMAKE_BINARY_DIR}/cmake_info.f90)
