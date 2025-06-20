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
    INSTALL_COMMAND make install && rm -rf ${CMAKE_BINARY_DIR}/external/greenx/src/greenx/.git/
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
)
set(GREENX_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/include/greenX/GNU-10.2.1/common/modules)
set(GREENX_LIBRARY_DIR ${CMAKE_BINARY_DIR}/external/greenx/install/lib)
include_directories(${GREENX_INCLUDE_DIR})
link_directories(${GREENX_LIBRARY_DIR})
