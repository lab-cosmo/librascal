/**
 * @file   bind_py_math.cc
 *
 * @author Felix Musil <musil.felix@gmail.com>
 *
 * @date   22 August 2018
 *
 * @brief  File for binding utils subroutines
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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


#include "bind_include.hh"



void math_binding(py::module& m) {
    m.def("hyp2f1", &math::hyp2f1, "y = hyp2f1( a, b, c, x )");
}
