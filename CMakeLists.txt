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

# external dependencies
include(ExternalProject)
include(cmake/fftw_dependency.cmake)
include(cmake/greenx_dependency.cmake)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(LIBS fftw3 gx_ac ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Compiler flags
set(CMAKE_Fortran_FLAGS "-Og -g -fopenmp -ffree-line-length-none -fcheck=bounds,do,recursion,pointer -cpp -I${CMAKE_BINARY_DIR}")

# Source files
set(SOURCES
    src/kinds.f90
    src/config_info.f90
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

# Website 
option(ENABLE_DOCS "enable website generation" OFF)
message(STATUS "Enable website generation: ${ENABLE_DOCS}")
if (${ENABLE_DOCS})
    include(cmake/website.cmake)
endif()

# save cmake variables in f90 file to make them accesseble from fortran
configure_file(src/cmake_info.f90.in ${CMAKE_BINARY_DIR}/cmake_info.f90)