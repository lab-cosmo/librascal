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

using namespace rascal; // NOLINT

template<typename RepresentationManager>
decltype(auto) add_representation_manager(py::module & mod, py::module & ) {
  using Manager_t = typename RepresentationManager::Manager_t;

  std::string representation_name =
        internal::GetBindingTypeName<RepresentationManager>();

  py::class_<RepresentationManager,
             RepresentationManagerBase>
             representation(mod, representation_name.c_str());
  // use custom constructor to pass json formated string as initializer
  // an alternative would be to convert python dict to json internally
  // but needs some workon in the pybind machinery
  representation.def(py::init([](Manager_t & manager, std::string& hyper_str) {
        // convert to json
        json hypers = json::parse(hyper_str);
        return std::make_unique<RepresentationManager>(manager, hypers);
                        }));
  representation.def("compute", &RepresentationManager::compute);

  return representation;
}


//! Representation python binding
void add_representation_managers(py::module & mod, py::module & m_garbage) {
  py::class_<RepresentationManagerBase>(m_garbage, "RepresentationManagerBase");
  using Manager_t = AdaptorStrict<AdaptorNeighbourList<
                                                StructureManagerCenters>>;
  using Representation1_t =
        RepresentationManagerSortedCoulomb<Manager_t>;
  auto rep_sorted_coulomb = add_representation_manager<
                                    Representation1_t>(mod, m_garbage);
  using Representation3_t =
          RepresentationManagerSphericalExpansion<Manager_t>;
  auto rep_spherical_expansion = add_representation_manager<Representation3_t>(
    mod, m_garbage);
}
