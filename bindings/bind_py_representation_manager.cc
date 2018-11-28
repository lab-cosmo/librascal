/**
 * @file   bind_py_representation_manager.cc
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
decltype(auto) add_representation_manager(py::module & mod, py::module & ) {
  using Manager_t = typename RepresentationManager::Manager_t;

  std::string representation_name =
        internal::GetBindingTypeName<RepresentationManager>();

  py::class_<RepresentationManager,
             RepresentationManagerBase>
             representation(mod, representation_name.c_str());
  representation.def(py::init<Manager_t &, std::string  >());
  representation.def("compute", &RepresentationManager::compute);

  return representation;
}


//! Representation python binding
void add_representation_managers(py::module & mod, py::module & m_garbage) {
  py::class_<RepresentationManagerBase>(m_garbage, "RepresentationManagerBase");
  using Manager_t = AdaptorStrict<AdaptorNeighbourList<
                                                StructureManagerCenters>>;
  using Representation_t =
        RepresentationManagerSortedCoulomb<Manager_t, CMoptions::Distance>;

  auto rep_sorted_coulomb = add_representation_manager<
                                    Representation_t>(mod, m_garbage);
}
