/**
 * @file   bind_py_neighbour_manager.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   9 Mai 2018
 *
 * @brief  File for binding the Neighbour Managers
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





#include <pybind11/pybind11.h>
#include "neighbourhood_managers/neighbourhood_manager_cell.hh"
#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <basic_types.hh>

using namespace rascal;
namespace py=pybind11;


void add_manager_cell(py::module& m){
	m.doc()  = "binding for the Neighbourhood Manager Linked Cell" ;
	py::class_<NeighbourhoodManagerBase<NeighbourhoodManagerCell>>(m, "NeighbourhoodManagerBase_Cell")
    .def(py::init<>());

	py::class_<NeighbourhoodManagerCell,NeighbourhoodManagerBase<NeighbourhoodManagerCell>>(m, "NeighbourhoodManagerCell")
    .def(py::init<>())
		.def("update",&NeighbourhoodManagerCell::update)
		//.def("build",[]())
		.def("__iter__", 
			[](NeighbourhoodManagerBase<NeighbourhoodManagerCell> &v) {
       	return py::make_iterator(v.begin(),v.end());
    	}, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */


	//py::class_<NeighbourhoodManagerBase<NeighbourhoodManagerCell>::iterator<1, 2>>(m, "Center_iterator")
	//	.def(py::init<>())

	py::class_<NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>>(m, "Cell.Center")
		.def("get_atom_index",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>::get_atom_index)
		.def("get_atom_type",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>::get_atom_type)
		.def("get_index",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>::get_index) 
		.def("get_size",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>::size)
		.def("get_position",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2>::get_position)
		.def("__iter__", 
			[](NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<1, 2> &v) {
       	return py::make_iterator(v.begin(),v.end());
    	}, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */

	py::class_<NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>>(m, "Cell.Neighbour")
		.def("get_atom_index",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::get_atom_index)
		.def("get_atom_type",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::get_atom_type)
		.def("get_index",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::get_index) 
		.def("get_size",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::size)
		//.def("get_atom_shift",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::get_atom_shift)
		.def("get_position",&NeighbourhoodManagerBase<NeighbourhoodManagerCell>::ClusterRef<2, 2>::get_position); 
}
