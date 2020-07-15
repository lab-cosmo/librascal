/**
 * @file   bind_py_math.hh
 *
 * @author
 *
 * @date   2020
 *
 * @brief  File for binding utils subroutines
 *
 * Copyright  , COSMO (EPFL), LAMMM (EPFL)
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

#include <pybind11/numpy.h>

#include "bind_include.hh"

#include "rascal/math/spherical_harmonics.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"

namespace rascal {
  void add_math(py::module &);
}

#endif  // BINDINGS_BIND_PY_MATH_HH_
