/**
 * @file   bind_py_math.hh
 *
 * @author Felix Musil <musil.felix@gmail.com>
 *
 * @date   22 August 2018
 *
 * @brief  File for binding utils subroutines
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef BINDINGS_BIND_PY_MATH_HH_
#define BINDINGS_BIND_PY_MATH_HH_

#include "math/math_interface.hh"
#include "math/math_utils.hh"

#include "bind_include.hh"

namespace rascal {
  void math_binding(py::module &);
}

#endif  // BINDINGS_BIND_PY_MATH_HH_
