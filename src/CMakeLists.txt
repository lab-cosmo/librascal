# =============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
# @author Felix Musil <felix.musil@epfl.ch>
# @author Markus Stricker <markus.stricker@epfl.ch>
#
# @date   01 Mar 2018
#
# @brief  configuration for system implementations
#
# @section LICENSE
#
# Copyright  2018 Till Junge, Felix Musil, Markus Stricker COSMO (EPFL),
# LAMMM (EPFL)
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

set(rascal_sources
    json_io.cc
    units.cc
    rascal_utility.cc
    math/math_utils.cc
    utils/sparsify_fps.cc
    structure_managers/structure_manager_lammps.cc
    structure_managers/structure_manager_centers.cc
    representations/calculator_base.cc
)

add_library(${LIBRASCAL_NAME} ${rascal_sources})

target_include_directories(${LIBRASCAL_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(${LIBRASCAL_NAME} PRIVATE -Werror)

target_link_libraries(${LIBRASCAL_NAME} PUBLIC Eigen3::Eigen)
target_link_libraries(${LIBRASCAL_NAME} PUBLIC ${WIGXJPF_NAME})

install(TARGETS ${LIBRASCAL_NAME} DESTINATION lib)
