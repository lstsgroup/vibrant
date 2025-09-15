# Check if GreenX paths are set in the initial cache file
if(DEFINED GREENX_INCLUDE_DIR AND DEFINED GREENX_LIBRARY_DIR)
    # Use the provided paths
    message(STATUS "Using provided GreenX paths:")
    message(STATUS "  Include directory: ${GREENX_INCLUDE_DIR}")
    message(STATUS "  Library directory: ${GREENX_LIBRARY_DIR}")
    include_directories(${GREENX_INCLUDE_DIR})
    link_directories(${GREENX_LIBRARY_DIR})

    # Define a dummy target for GreenX
    add_custom_target(greenx)
else()
    # Build GreenX as an external project
    message(STATUS "GreenX paths not provided. Building GreenX as an external project.")
    ExternalProject_Add(greenx
        GIT_REPOSITORY https://github.com/nomad-coe/greenX.git
        GIT_TAG v2.0  # You can specify a tag or branch here
        PREFIX ${CMAKE_BINARY_DIR}/external/greenx
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/greenx/install
                   -DMINIMAX_COMPONENT=OFF
                   -DLBASIS_COMPONENT=OFF
                   -DPAW_COMPONENT=OFF
                   -DCOMPILE_SUBMODULES=OFF
        BUILD_COMMAND make -j 
        INSTALL_COMMAND make install && rm -rf ${CMAKE_BINARY_DIR}/external/greenx/src/greenx/.git/
        GIT_SUBMODULES ""
        UPDATE_COMMAND ""
    )
    set(GREENX_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/include/modules)
    set(GREENX_LIBRARY_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/lib)
    message(STATUS "GreenX include directory set to: ${GREENX_INCLUDE_DIR}")
    message(STATUS "GreenX library directory set to: ${GREENX_LIBRARY_DIR}")
    include_directories(${GREENX_INCLUDE_DIR})
    link_directories(${GREENX_LIBRARY_DIR})
endif()
