# python module for building the website
find_package(Sphinx REQUIRED)

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

#add_custom_target(docs ALL
#    ${SPHINX_EXECUTABLE}
#        -b html
#        -d "${SPHINX_CACHE_DIR}"
#        "${CMAKE_CURRENT_SOURCE_DIR}/docs/"
#        "${SPHINX_HTML_DIR}"
#    COMMENT "Building HTML documentation with Sphinx")
add_custom_target(docs ALL
    COMMAND ${SPHINX_EXECUTABLE}
        -b html
        -c "${CMAKE_CURRENT_SOURCE_DIR}/docs"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/docs"
        "${SPHINX_HTML_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs"
    COMMENT "Building HTML documentation with Sphinx")