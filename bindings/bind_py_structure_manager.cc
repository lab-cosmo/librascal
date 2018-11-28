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


#include "structure_managers/structure_manager.hh"
#include "structure_managers/structure_manager_centers.hh"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include <basic_types.hh>

using namespace rascal;
namespace py = pybind11;

template<size_t Order, typename StructureManagerImplementation>
decltype(auto) add_cluster(py::module & m) {
  using ClusterRef = typename StructureManager<
    StructureManagerImplementation>::template ClusterRef<Order>;
  // TODO: change the exposed name to a convertion of the
  // StructureManagerImplementation type to a string
  py::class_<ClusterRef>
    py_cluster(m, (Order == 1)
                ? "StructureManager.Center" : "StructureManager.Neighbour");
  py_cluster.def_property_readonly("atom_index", & ClusterRef::get_atom_index,
                                   py::return_value_policy::reference)
    .def_property_readonly("atom_type", & ClusterRef::get_atom_type,
                           py::return_value_policy::reference)
    .def_property_readonly("index", & ClusterRef::get_index,
                           py::return_value_policy::reference)
    .def_property_readonly("size", & ClusterRef::size,
                           py::return_value_policy::reference)
    .def_property_readonly("position", & ClusterRef::get_position,
                           py::return_value_policy::reference);
  return py_cluster;
}

void add_manager_centers(py::module & m){
  m.doc() = "binding for the Structure Manager Centers";
  py::class_<StructureManager<StructureManagerCenters>>
    (m, "StructureManagerBase_Centers").def(py::init<>());
  py::class_<StructureManagerCenters,
             StructureManager<StructureManagerCenters>>
    (m, "StructureManagerCenters")
    .def(py::init<>())
    //! interface can handle both row- and col-major matrices
    .def("update", [] (StructureManagerCenters & v,
                       const py::EigenDRef<const Eigen::MatrixXd> & positions,
                       const py::EigenDRef<const Eigen::VectorXi> & atom_types,
                       const py::EigenDRef<const Eigen::MatrixXd> & cell,
                       const py::EigenDRef<const Eigen::MatrixXi> & pbc ) {
           v.update(positions, atom_types, cell, pbc);
         })
    .def("__iter__", [] (StructureManager<StructureManagerCenters> & v) {
        return py::make_iterator(v.begin(), v.end());
      }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
  add_cluster<1, StructureManagerCenters>(m);
}
