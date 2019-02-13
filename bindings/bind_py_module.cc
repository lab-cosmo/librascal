/**
 * @file   bind_py_module.cc
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  Main binding file for Rascal
 *
 * Copyright © 2017 Felix Musil
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

using namespace pybind11::literals;  // NOLINT (is recommended use of pybind11)
namespace py = pybind11;

PYBIND11_MODULE(_rascal, mod) {
  mod.doc() = "Python bindings for the Rascal library";

  py::module m_nl = mod.def_submodule("NeighbourList");
  m_nl.doc() = "Utilities to make neighbour lists";

  py::module m_rpr_mng = mod.def_submodule("RepresentationManager");
  m_rpr_mng.doc() = "Representation Manager Classes";
  py::module m_feat_mng = mod.def_submodule("FeatureManager");
  m_feat_mng.doc() = "Feature Manager Classes";

  py::module m_utl = mod.def_submodule("utils");
  py::module m_math = mod.def_submodule("math");
  m_math.doc() = "Collection of math functions";
  py::module m_garbage = mod.def_submodule("rubbish");
  m_garbage.doc() = "Collection of bindings that are needed but not functional";

  py::add_ostream_redirect(m_utl, "ostream_redirect");

  add_structure_managers(m_nl, m_garbage);
  add_representation_managers(m_rpr_mng, m_garbage);
  add_feature_managers(m_feat_mng, m_garbage);
  utils_binding(m_utl);
  math_binding(m_math);
}
