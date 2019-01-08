# =============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
#
# @date   08 Jan 2018
#
# @brief  configuration for system implementations
#
# @section LICENSE
#
# Copyright © 2018 Till Junge
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
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_representation_manager.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_feature_manager.cc
  )



install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/rascal/" DESTINATION ./
          PATTERN "*__pycache__" EXCLUDE
          PATTERN "*.so" EXCLUDE
          PATTERN "*.dylib" EXCLUDE)
# see for more details on the linking path handling in Linux and MacOS when
# installing:
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
# forces the writing of absolute rpath into the library when installing
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON")
# set(CMAKE_SKIP_RPATH "OFF")
#set(CMAKE_MACOSX_RPATH "ON")
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/")

# builds the python bindings (must be after setting up the rpath)
pybind11_add_module(pyRascal ${PY_BINDING_SRCS})
target_link_libraries(pyRascal PRIVATE rascal cephes)
# Want to rename the output, so that the python module is called _rascal
set_target_properties(pyRascal PROPERTIES OUTPUT_NAME _rascal)
# Sets rule to install the bindings library in the rascal folder
install(TARGETS pyRascal DESTINATION lib)

# Copy rascal in the build folder as part of the make process of bindings build
add_custom_command(TARGET pyRascal POST_BUILD
  COMMAND cp -R
  ${CMAKE_CURRENT_SOURCE_DIR}/rascal ${CMAKE_BINARY_DIR}
  )
# Copy the resulting bindings lib in the python wrapper lib folder
add_custom_command(TARGET pyRascal POST_BUILD
  COMMAND cp
  ${CMAKE_BINARY_DIR}/bindings/_rascal.* ${CMAKE_BINARY_DIR}/rascal/lib/
  )
