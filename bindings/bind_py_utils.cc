/**
 * @file   bind_py_utils.cc
 *
 * @author Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   22 August 2018
 *
 * @brief  File for binding utils subroutines
 *
 * Copyright  2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
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

#include "bind_py_utils.hh"

namespace rascal {
  void utils_binding(py::module & mod) {
    py::module m_utils_sparse = mod.def_submodule("sparsification");
    m_utils_sparse.doc() = "Sparsification Routines";

    m_utils_sparse.def(
        "fps", &utils::select_fps,
        "Selects points from a NxD dimensional"
        " feature matrix by farthest point sampling (N is the number of"
        " sample in a D dimensional space).",
        py::arg("feature_matrix"), py::arg("n_sparse") = 0,
        py::arg("i_first_point") = 0,
        py::arg("restart") = std::make_tuple(
            Eigen::ArrayXi(0), Eigen::ArrayXd(0), Eigen::ArrayXd(0)));

    m_utils_sparse.def(
        "fps_voronoi", &utils::select_fps_voronoi,
        "Selects points from a"
        " NxD dimensional feature matrix by farthest point sampling, using"
        " a Voronoi cell method (N is the number of sample in a D dimensional"
        " space).",
        py::arg("feature_matrix"), py::arg("n_sparse"),
        py::arg("i_first_point"));
  }
}  // namespace rascal
