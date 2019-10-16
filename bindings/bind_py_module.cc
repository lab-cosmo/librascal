/**
 * @file   bind_py_module.cc
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  Main binding file for Rascal
 *
 * Copyright  2017 Felix Musil
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

#include "bind_py_module.hh"

using namespace pybind11::literals;  // NOLINT (is recommended use of pybind11)
namespace py = pybind11;

PYBIND11_MODULE(_rascal, mod) {
  mod.doc() = "Python bindings for the Rascal library";

  py::module m_nl = mod.def_submodule("neighbour_list");
  m_nl.doc() = "Utilities to make neighbour lists";

  py::module m_repr = mod.def_submodule("representation_calculators");
  m_repr.doc() = "Representation calculator classes";
  py::module m_models = mod.def_submodule("models");
  m_models.doc() = "Collection of models";

  py::module m_kernels = m_models.def_submodule("kernels");
  m_kernels.doc() = "Collection of Kernels";

  py::module m_utl = mod.def_submodule("utils");
  py::module m_internal = mod.def_submodule("_internal");
  m_internal.doc() =
      "Collection of bindings that are needed to build functional bindings for"
      " the python user but are not functional itself. It basically contains"
      " all methods which are used on the binding side, but are not meant to"
      " be used by the python user. It is also not part of the rascal library";

  py::add_ostream_redirect(m_utl, "ostream_redirect");

  rascal::add_structure_managers(m_nl, m_internal);
  rascal::add_representation_calculators(m_repr, m_internal);
  rascal::utils_binding(m_utl);
  rascal::add_kernels(m_kernels, m_internal);
}
