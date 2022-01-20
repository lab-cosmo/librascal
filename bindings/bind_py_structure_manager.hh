/**
 * @file   bind_py_structure_manager.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   9 Mai 2018
 *
 * @brief  File for binding the Neighbour Managers
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef BINDINGS_BIND_PY_STRUCTURE_MANAGER_HH_
#define BINDINGS_BIND_PY_STRUCTURE_MANAGER_HH_

#include "bind_include.hh"

#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_full_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_half_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_increase_maxorder.hh"
#include "rascal/structure_managers/adaptor_kspace.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager.hh"
#include "rascal/structure_managers/structure_manager_base.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/structure_managers/structure_manager_lammps.hh"

namespace rascal {
  void add_structure_managers(py::module &, py::module & /*m_unused*/);
}
#endif  // BINDINGS_BIND_PY_STRUCTURE_MANAGER_HH_
