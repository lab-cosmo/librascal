/**
 * @file   bind_py_utils.cc
 *
 * @author Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   22 August 2018
 *
 * @brief  File for binding utils subroutines
 *
 * Copyright © 2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
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


#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "utils/sparsify_utilities.hh"

using namespace rascal;
namespace py=pybind11;


void utils_binding(py::module& m){
    m.def("fps", &utils::fps, "Selects points from a NxD dimensional feature matrix by farthest point sampling",
          py::arg("feature_matrix"), py::arg("n_sparse"), py::arg("i_first_point") );
}
