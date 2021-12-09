/**
 * @file   bind_py_representation_calculator.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   30 Oct 2018
 *
 * @brief  File for binding the Representation calculator
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_
#define BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_

#include "bind_include.hh"
#include "bind_py_structure_manager.hh"

#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_covariants.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"

namespace rascal {
  void add_representation_calculators(py::module &, py::module &);
}

#endif  // BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_
