find_package(Sphinx REQUIRED)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/Doxyfile"
  "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/input_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/input_def"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/source/xml_output_def.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xml_output_def"
  @ONLY)

add_custom_target(dev_doc
  ${SPHINX_EXECUTABLE}
  -q -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/source"
  "${SPHINX_HTML_DIR}"
  COMMENT "Building HTML documentation with Sphinx")

