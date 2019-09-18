/**
 * @file   models.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   16 Jun 2019
 *
 * @brief  File for binding the models
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef BINDINGS_BIND_PY_MODELS_HH_
#define BINDINGS_BIND_PY_MODELS_HH_

#include "bind_include.hh"
#include "bind_py_structure_manager.hh"
#include "bind_py_representation_calculator.hh"

#include "models/kernels.hh"

namespace rascal {
  void add_kernels(py::module &, py::module &);
}

#endif  // BINDINGS_BIND_PY_MODELS_HH_
