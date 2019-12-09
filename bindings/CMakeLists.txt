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

set(BINDING_SRCS
    bind_py_module.cc
    bind_py_structure_manager.cc
    bind_py_utils.cc
    bind_py_representation_calculator.cc
    bind_py_models.cc
)

# copy rascal folder in the python distribution
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/rascal" DESTINATION "${PYTHON_SITE_PACKAGES}"
    PATTERN "*__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
    PATTERN "*.so" EXCLUDE
    PATTERN "*.dylib" EXCLUDE
)

# builds the python bindings (must be after setting up the rpath)
# SYSTEM: use include_directory(SYSTEM)
# NO_EXTRAS: do not setup LTO (link time optimization) and striping of the
#            shared library 
pybind11_add_module(_rascal ${BINDING_SRCS} SYSTEM NO_EXTRAS)
target_compile_options(_rascal PRIVATE -Werror)

target_link_libraries(_rascal PRIVATE "${LIBRASCAL_NAME}")
# Remove the 'lib' prefix, so that the python module is called _rascal
set_target_properties(_rascal PROPERTIES PREFIX "")
# Sets rule to install the bindings library in the rascal folder
install(TARGETS _rascal DESTINATION "${PYTHON_SITE_PACKAGES}/rascal/lib")

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
