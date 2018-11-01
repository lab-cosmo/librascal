/**
 * @file   bind_include.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   31 Oct 2018
 *
 * @brief  File to centralize includes and function declaration
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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



#ifndef BIND_INCLUDE_H
#define BIND_INCLUDE_H

#include "utils/sparsify_utilities.hh"

#include "math/math_interface.hh"
#include "math/math_utils.hh"

#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_sorted_coulomb.hh"

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/structure_manager_lammps.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/structure_manager.hh"

#include "basic_types.hh"
#include "rascal_utility.hh"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include <Eigen/Dense>
#include <vector>


// PYBIND11_MAKE_OPAQUE(std::vector<double>);

using namespace rascal;
namespace py=pybind11;

void add_structure_managers(py::module&,py::module&,py::module&);
void add_representation_managers(py::module&,py::module&);

void utils_binding(py::module&);
void math_binding(py::module&);



#endif /* BIND_INCLUDE_H */