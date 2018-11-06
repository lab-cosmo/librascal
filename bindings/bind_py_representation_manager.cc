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




template<typename RepresentationManager>
decltype(auto) add_representation_manager(py::module & mod,py::module & ){
  using Manager_t = typename RepresentationManager::Manager_t;

  std::string representation_name = 
        internal::GetBindingTypeName<RepresentationManager>();

  py::class_<RepresentationManager, 
             RepresentationManagerBase> 
             representation(mod, representation_name.c_str());
  representation.def(py::init<Manager_t &, double , double , double , size_t  >());
  representation.def("compute", &RepresentationManager::compute);
  
  return representation;
};


//! Sorted Coulomb representation python binding
void add_representation_managers(py::module & mod,py::module & m_garbage){
  
  py::class_<RepresentationManagerBase>(m_garbage,"RepresentationManagerBase")
      .def(py::init<>());

  using Representation_t = RepresentationManagerSortedCoulomb<
        AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;
  
  auto representation = add_representation_manager<
                                    Representation_t>(mod,m_garbage);
  // TODO will have to change to an external data structure handler
  representation.def("get_representation_full", 
         &Representation_t::get_representation_full,
         py::return_value_policy::reference_internal,py::keep_alive<1,0>());
  representation.def("get_coulomb_matrix",
          &Representation_t::get_coulomb_matrix<1,1>,
          py::return_value_policy::reference_internal,py::keep_alive<1,0>());
};
