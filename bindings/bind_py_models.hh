/**
 * @file   bind_py_models.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   16 Jun 2019
 *
 * @brief  File for binding the models
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef BINDINGS_BIND_PY_MODELS_HH_
#define BINDINGS_BIND_PY_MODELS_HH_

#include "bind_include.hh"
#include "bind_py_representation_calculator.hh"
#include "bind_py_structure_manager.hh"

#include "rascal/models/kernels.hh"
#include "rascal/models/numerical_kernel_gradients.hh"
#include "rascal/models/sparse_kernel_predict.hh"
#include "rascal/models/sparse_kernels.hh"
#include "rascal/models/sparse_points.hh"

namespace rascal {
  void add_models(py::module &, py::module &);
}

#endif  // BINDINGS_BIND_PY_MODELS_HH_
