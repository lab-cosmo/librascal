/**
 * @file   bind_py_cdist.cc
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  File for binding the function cdist to a python object
 *
 * Copyright © 2017 Felix Musil
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
#include "cdist.hh"
#include <pybind11/eigen.h>

using namespace rascal;
namespace py=pybind11;

void add_cdist(py::module& m)
{
	m.doc()  = " binding for the distance matrix calculation" ;
	m.def("cdist",&cdist);
	//m.def("scale", [](py::EigenDRef<Eigen::MatrixXd> mm, double c) { mm *= c; });
}

