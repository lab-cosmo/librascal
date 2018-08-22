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
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# Rascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
# =============================================================================



set (PY_BINDING_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_module.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_structure_manager.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bind_py_utils.cc
  )

find_package(PythonLibsNew 3 MODULE REQUIRED)

pybind11_add_module(pyRascal ${PY_BINDING_SRCS})
target_link_libraries(pyRascal PRIVATE rascal)
# Want to rename the output, so that the python module is called _rascal
set_target_properties(pyRascal PROPERTIES OUTPUT_NAME _rascal)
