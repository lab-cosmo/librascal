find_package(Sphinx REQUIRED)
find_package(Doxygen REQUIRED)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

set(DOXYGEN_BUILD_HTML "NO" CACHE BOOL
    "Build the native Doxygen HTML documentation \
(in addition to sphinx/breathe)?")

mark_as_advanced(DOXYGEN_BUILD_HTML)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/Doxyfile"
  "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
  @ONLY)

# Copy over the resources needed for the Doxygen HTML docs (transitional)
file(GLOB DOX_RESOURCE_FILES
          "${CMAKE_CURRENT_SOURCE_DIR}/source/dox_resources/*")
file(GLOB DOX_MARKDOWN_FILES
          "${CMAKE_CURRENT_SOURCE_DIR}/source/*.md")
file(COPY ${DOX_RESOURCE_FILES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY ${DOX_MARKDOWN_FILES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/input_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/input_def"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/xml_output_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xml_output_def"
  @ONLY)

add_custom_target(doxygen_doc
    COMMAND "${DOXYGEN_EXECUTABLE}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating C++ documentation with Doxygen")

add_custom_target(doc
  COMMAND
  ${SPHINX_EXECUTABLE}
  -q -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/source"
  "${SPHINX_HTML_DIR}"
  DEPENDS doxygen_doc
  COMMENT "Building HTML documentation with Sphinx")
