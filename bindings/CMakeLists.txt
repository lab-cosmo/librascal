# =============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
# @author Felix Musil <felix.musil@epfl.ch>
#
# @date   08 Jan 2018
#
# @brief  configuration for system implementations
#
# @section LICENSE
#
# Copyright  2018 Till Junge
#
# Rascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# Rascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
# =============================================================================

set (PY_BINDING_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_module.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_structure_manager.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_utils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_math.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_representation_calculator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_models.cc
  )


# copy rascal folder in the python distribution
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/rascal" DESTINATION "${PYTHON_SITE_PACKAGES}"
          PATTERN "*__pycache__" EXCLUDE
          PATTERN "*.pyc" EXCLUDE
          PATTERN "*.so" EXCLUDE
          PATTERN "*.dylib" EXCLUDE)


# builds the python bindings (must be after setting up the rpath)
pybind11_add_module(pyRascal ${PY_BINDING_SRCS})
target_link_libraries(pyRascal PRIVATE "${LIBRASCAL_NAME}")
# Want to rename the output, so that the python module is called _rascal
set_target_properties(pyRascal PROPERTIES OUTPUT_NAME _rascal)
# Sets rule to install the bindings library in the rascal folder
install(TARGETS pyRascal DESTINATION "${PYTHON_SITE_PACKAGES}/rascal/lib")

# force cmake to run these commands at every make
# https://stackoverflow.com/questions/13920072/how-to-always-run-command-when-building-regardless-of-any-dependency
add_custom_target(librascal ALL
    DEPENDS
        copy_py
    )

add_dependencies(librascal pyRascal)

# Copy the python wrappers in the build folder
# Copy the resulting bindings lib in the python wrapper lib folder
add_custom_command(OUTPUT copy_py
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/rascal ${CMAKE_BINARY_DIR}/rascal
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_BINARY_DIR}/bindings/_rascal.* ${CMAKE_BINARY_DIR}/rascal/lib/
  DEPENDS
      always_rebuild
  )


add_custom_command(
    OUTPUT always_rebuild
    COMMAND ${CMAKE_COMMAND} -E echo
)

