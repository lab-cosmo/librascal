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
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_
#define BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_

#include "bind_include.hh"
#include "bind_py_structure_manager.hh"
#include "representations/calculator_base.hh"
#include "representations/calculator_sorted_coulomb.hh"
#include "representations/calculator_spherical_covariants.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"

namespace rascal {
  void add_representation_calculators(py::module &, py::module &);
}

#endif  // BINDINGS_BIND_PY_REPRESENTATION_CALCULATOR_HH_
