#==============================================================================
# file   CMakeLists.txt
#
# @author Felix Musil <felix.musil@epfl.ch>
#
# @date  14 Oct 2018
#
# @brief  Configuration for compiling all representations
#
# @section LICENSE
#
# Copyright © 2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
#
# librascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# librascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#

# locate all cephes subfolders to process
FILE(GLOB  cephes_SRC
    "bessel/*.c"
    "cmath/*.c"
    "cprob/*.c"
    "misc/*.c"
    "ellf/*.c"
    "polyn/*.c"
    "mtherr.c"
    )
# note: single/ is not compiled because it conflicts
# with many functions in cmath and cprob :(

# use our custom mtherr instead of cmath/mtherr.c
LIST(REMOVE_ITEM  cephes_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/cmath/mtherr.c"
    )


# ignore files which have a main()
LIST(REMOVE_ITEM  cephes_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/cmath/mtst.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmath/dtestvec.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmath/mod2pi.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ellf/ellf.c")

# ignore files who are exact duplicates of cmath
LIST(REMOVE_ITEM  cephes_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/cprob/drand.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cprob/const.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cprob/mtherr.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cprob/polevl.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/cprob/unity.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/chbevl.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/mtherr.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/polevl.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ellf/cmplx.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ellf/const.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ellf/mtherr.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ellf/polevl.c"
    )

# ignore files in misc who are exact duplicates of any other
LIST(REMOVE_ITEM  cephes_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/psi.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/revers.c"
    )

# ignore rational polynomials of polyr.c,
# because it conflict with the simpler polynomials of polyn.c
LIST(REMOVE_ITEM cephes_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/polyn/polyr.c"
    )

target_sources("${LIBCEPHES_NAME}" PRIVATE ${cephes_SRC})

