# Check if BLAS and LAPACK paths are provided
if(DEFINED BLAS_LIBRARY_DIR AND DEFINED LAPACK_LIBRARY_DIR)
    message(STATUS "Using provided BLAS/LAPACK paths:")
    message(STATUS "  BLAS library directory: ${BLAS_LIBRARY_DIR}")
    message(STATUS "  LAPACK library directory: ${LAPACK_LIBRARY_DIR}")
    link_directories(${BLAS_LIBRARY_DIR} ${LAPACK_LIBRARY_DIR})
    set(BLAS_LIBRARIES "${BLAS_LIBRARY_DIR}/libblas.so")  # Adjust the library name as needed
    set(LAPACK_LIBRARIES "${LAPACK_LIBRARY_DIR}/liblapack.so")  # Adjust the library name as needed
else()
    # Use find_package to locate BLAS and LAPACK
    message(STATUS "BLAS/LAPACK paths not provided. Using find_package to locate them.")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()
