cmake_minimum_required(VERSION 3.18)

# Project name and Fortran standard
project(vibrant Fortran)

include(ExternalProject)

# Build FFTW as an external project
ExternalProject_Add(fftw
    URL https://www.fftw.org/fftw-3.3.10.tar.gz
    PREFIX ${CMAKE_BINARY_DIR}/external/fftw
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/external/fftw/install
    BUILD_COMMAND make -j 
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
)
set(FFTW_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/fftw/install/include)
set(FFTW_LIBRARY_DIR ${CMAKE_BINARY_DIR}/external/fftw/install/lib)

# Build GreenX as an external project
ExternalProject_Add(greenx
    GIT_REPOSITORY https://github.com/nomad-coe/greenX.git
    GIT_TAG main  # You can specify a tag or branch here
    PREFIX ${CMAKE_BINARY_DIR}/external/greenx
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/greenx/install
               -DMINIMAX_COMPONENT=OFF
               -DLBASIS_COMPONENT=OFF
               -DPAW_COMPONENT=OFF
               -DCOMPILE_SUBMODULES=OFF
    BUILD_COMMAND make -j 
    INSTALL_COMMAND make install
    GIT_SUBMODULES ""
)
set(GREENX_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/include/modules)
set(GREENX_LIBRARY_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/lib)

# Compiler flags
set(CMAKE_Fortran_FLAGS "-Og -g -fopenmp -fcheck=bounds,do,recursion,pointer")
include_directories(${FFTW_INCLUDE_DIR} ${GREENX_INCLUDE_DIR})

# Linker flags
link_directories(${FFTW_LIBRARY_DIR} ${GREENX_LIBRARY_DIR})

# Libraries to link
set(LIBS fftw3 gx_ac)

# Source files
set(SOURCES
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


option(ENABLE_REGRESSION_TESTS "Enable vibrant regreesion tests" On)
if (${ENABLE_REGRESSION_TESTS})
    # Python required for application testing
    include(cmake/python3.cmake)
    find_package(Python3 3.10 COMPONENTS Interpreter Development)
    if (Python3_FOUND)
        message("-- Python 3 interpreter version: " ${Python3_VERSION})
    else()
        message("-- Python 3 interpreter not found")
    endif()

    # Enable ctest
    enable_testing()

    # copy the test to build directory
    add_custom_command(
	    TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/test
                ${PROJECT_BINARY_DIR}/test)

    add_test(
        NAME IR_Berry
        COMMAND pytest test_IR_Berry.py
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/IR_Berry
    )
endif()

