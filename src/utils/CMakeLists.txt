#==============================================================================
# file   CMakeLists.txt
#
# @author Michele Ceriotti <michele.ceriotti@gmail.com>
#
# @date   22 August 2018
#
# @brief  configuration for compiling utilities
#
# @section LICENSE
#
# Copyright © 2018 Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
#
# rascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# rascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#


set (utils_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/sparsify_fps.cc
  )

target_sources(rascal PRIVATE ${utils_SRC})
