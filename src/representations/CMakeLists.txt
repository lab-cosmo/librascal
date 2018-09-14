#==============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
#
# @date   04 May 2018
#
# @brief  Configuration for compiling all representations
#
# @section LICENSE
#
# Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
#
# librascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# librascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#


set(representation_SRC
  representations/representation_manager_base.cc
  )

target_sources(rascal PRIVATE ${representation_SRC})
