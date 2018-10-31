/**
 * @file   bind_py_neighbour_manager.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   30 Oct 2018
 *
 * @brief  File for binding the Representation Managers
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



#include "bind_include.hh"


//! Sorted Coulomb representation python binding
void add_sorted_coulomb(py::module & m){
  m.doc() = "binding for the Sorted Coulomb Representation" ;
  py::class_<RepresentationManagerBase>(m,"RepresentationManagerBase")
      .def(py::init<>());

  using Manager_t = AdaptorStrict<AdaptorNeighbourList<
                            StructureManagerCenters>>;
  
  py::class_<RepresentationManagerSortedCoulomb<Manager_t>,
              RepresentationManagerBase> (m, 
                      "SortedCoulombRepresentation")
    .def(py::init<Manager_t &, double , 
      double , double , 
      size_t  >())
    .def("compute",
        &RepresentationManagerSortedCoulomb<Manager_t>::compute)
    .def("get_representation_full",
        &RepresentationManagerSortedCoulomb<Manager_t>::get_representation_full,
        py::return_value_policy::copy);
    //py::keep_alive<1,0>(),
        
};
