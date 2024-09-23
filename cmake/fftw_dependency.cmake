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
include_directories(${FFTW_INCLUDE_DIR})
link_directories(${FFTW_LIBRARY_DIR})