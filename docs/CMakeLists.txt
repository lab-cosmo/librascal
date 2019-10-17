if (NOT ${BUILD_BINDINGS})
    message(FATAL_ERROR "BUILD_BINDINGS=ON is required to build the documentation")
endif()

find_package(Sphinx REQUIRED)
find_package(Doxygen REQUIRED)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

option(DOXYGEN_BUILD_HTML "Build the native Doxygen HTML documentation" OFF)
mark_as_advanced(DOXYGEN_BUILD_HTML)

if(${DOXYGEN_BUILD_HTML})
    set(DOXYGEN_GENERATE_HTML "YES")
else()
    set(DOXYGEN_GENERATE_HTML "NO")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
  @ONLY)

file(GLOB_RECURSE RASCAL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../src/**.hh)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygenxml/index.xml
    COMMAND "${DOXYGEN_EXECUTABLE}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating C++ documentation with Doxygen"
    DEPENDS ${RASCAL_HEADERS}
)

add_custom_target(doxygen_doc
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygenxml/index.xml
)

add_custom_target(doc ALL
  COMMAND ${SPHINX_EXECUTABLE} -q -b html
            -c "${CMAKE_CURRENT_BINARY_DIR}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}/source"
            "${SPHINX_HTML_DIR}"
  DEPENDS doxygen_doc copy_python_bindings
  COMMENT "Building HTML documentation with Sphinx"
)
