# Check if FFTW paths are set in the initial cache file
if(DEFINED FFTW_INCLUDE_DIR AND DEFINED FFTW_LIBRARY_DIR)
    # Use the provided paths
    message(STATUS "Using provided FFTW paths:")
    message(STATUS "  Include directory: ${FFTW_INCLUDE_DIR}")
    message(STATUS "  Library directory: ${FFTW_LIBRARY_DIR}")
    include_directories(${FFTW_INCLUDE_DIR})
    link_directories(${FFTW_LIBRARY_DIR})

    # Define a dummy target for FFTW
    add_custom_target(fftw)
else()
    # Build FFTW as an external project
    message(STATUS "FFTW paths not provided. Building FFTW as an external project.")
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
    message(STATUS "FFTW include directory set to: ${FFTW_INCLUDE_DIR}")
    message(STATUS "FFTW library directory set to: ${FFTW_LIBRARY_DIR}")
    include_directories(${FFTW_INCLUDE_DIR})
    link_directories(${FFTW_LIBRARY_DIR})
endif()