# =============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
#
# @date   01 Mar 2018
#
# @brief  configuration for system implementations
#
# @section LICENSE
#
# Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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



add_library(rascal "") # ${rascal_SRC})
target_link_libraries(rascal PUBLIC Eigen3::Eigen)

add_subdirectory(structure_managers)
add_subdirectory(representations)
add_subdirectory(utils)

install(TARGETS rascal DESTINATION lib)
