set(BINDING_SRCS
    bind_py_module.cc
    bind_py_structure_manager.cc
    bind_py_representation_calculator.cc
    bind_py_models.cc
)

# builds the python bindings (must be after setting up the rpath)
# SYSTEM: use include_directory(SYSTEM)
# NO_EXTRAS: do not setup LTO (link time optimization) and striping of the
#            shared library
pybind11_add_module(_rascal ${BINDING_SRCS} SYSTEM NO_EXTRAS)

target_link_libraries(_rascal PRIVATE "${LIBRASCAL_NAME}")
# Remove the 'lib' prefix, so that the python module is called _rascal
set_target_properties(_rascal PROPERTIES PREFIX "")
# Sets rule to install the bindings library in the rascal folder

# relative RPATH so that _rascal.so can find librascal.so when installed.
# see link for more details on the topic
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set_target_properties(_rascal PROPERTIES SKIP_BUILD_RPATH  FALSE)
set_target_properties(_rascal PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_target_properties(_rascal PROPERTIES BUILD_RPATH_USE_ORIGIN TRUE)
if(NOT APPLE AND UNIX)
  # this link path accommodates most installations
  # `$ORIGIN/../../../../` should allow to link libraries in `$PREFIX/lib`
  # `$ORIGIN/../../src` should allow to link libraries from the build folder
  # `$ORIGIN/../../lib` accommodates installation with
  # `python setup.py install'
  set(BINDNIGS_INSTALL_RPATH "$ORIGIN/../../../../" "$ORIGIN/../../lib"
                          "$ORIGIN/../../src" "$ORIGIN/../src")
elseif(APPLE)
  set(MACOSX_RPATH ON)
  # this link path accommodates most installations
  # `@loader_path/../../../../` should allow to link libraries
  # in `$PREFIX/lib`
  # `@loader_path/../../src` should allow to link libraries from the build folder
  # `@loader_path/../../lib` accommodates installation with
  # `python setup.py install'
  set(BINDNIGS_INSTALL_RPATH "@loader_path/../../../../" "@loader_path/../../lib"
                          "@loader_path/../../src" "@loader_path/../src")
endif()

set_target_properties(_rascal PROPERTIES INSTALL_RPATH "${BINDNIGS_INSTALL_RPATH}")

# Normal installation procedure
# copy rascal folder in the python distribution
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/rascal" DESTINATION "${CMAKE_INSTALL_PREFIX}/bindings"
    PATTERN "*__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
    PATTERN "*.so" EXCLUDE
    PATTERN "*.dylib" EXCLUDE
)
install(TARGETS _rascal DESTINATION "${CMAKE_INSTALL_PREFIX}/bindings/rascal/lib")

# Copy the python wrappers in the build folder
# Copy the resulting bindings lib in the python wrapper lib folder
add_custom_target(copy_python_bindings ALL
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/rascal
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/rascal ${CMAKE_BINARY_DIR}/rascal
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/_rascal.* ${CMAKE_BINARY_DIR}/rascal/lib/
    DEPENDS _rascal
    COMMENT "Copying python bindings and libraries"
)
