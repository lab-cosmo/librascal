#==============================================================================
# file   CMakeLists.txt
#
# @author Felix Musil <felix.musil@epfl.ch>
#
# @date   14 Oct 2018
#
# @brief  configuration for compiling interface to math library
#
# @section LICENSE
#
# Copyright © 2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

add_library(cephes "")

add_subdirectory(cephes)

install(TARGETS cephes DESTINATION lib)


add_library(rascal-math "")

install(TARGETS rascal-math DESTINATION lib)

set (math_SRC
  math_utils.cc
  math_interface.cc
  )

target_sources(rascal-math PRIVATE ${math_SRC})

target_link_libraries(rascal-math PRIVATE cephes)
