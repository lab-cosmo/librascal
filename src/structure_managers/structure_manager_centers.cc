/**
 * @file   structure_manager_centers.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   06 August 2018
 *
 * @brief Manager with atoms and centers
 *
 * Copyright  2018  Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "structure_managers/structure_manager_centers.hh"

#include <numeric>
#include <fstream>
#include <iostream>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  // function for setting the internal data structures
  void StructureManagerCenters::build() {
    auto && center_atoms_mask = this->get_center_atoms_mask();
    this->natoms = this->get_positions().size() / traits::Dim;
    this->n_center_atoms = center_atoms_mask.count();
    // initialize necessary data structure
    this->atoms_index[0].clear();
    this->offsets.clear();
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    // set the references to the center atoms positions and types
    for (size_t id{0}; id < this->natoms; ++id) {
      if (center_atoms_mask(id)) {
        this->atoms_index[0].push_back(id);
        this->offsets.push_back(id);
      }
    }

    for (size_t id{0}; id < this->natoms; ++id) {
      if (not center_atoms_mask(id)) {
        this->atoms_index[0].push_back(id);
        this->offsets.push_back(id);
      }
    }

    Cell_t lat = this->atoms_object.cell;
    this->lattice.set_cell(lat);

    // Check if all atoms are inside the unit cell assuming the cell starts
    // at (0,0,0)
    auto positions_scaled = this->atoms_object.get_scaled_positions();
    double tol{1e-10};
    if ((positions_scaled.array().rowwise().minCoeff() < -tol).any() or
        (positions_scaled.array().rowwise().maxCoeff() > 1. + tol).any()) {
      std::string error{R"(Some of the positions in the structure are not
                            inside the unit cell. Please wrap the atoms.)"};
      throw std::runtime_error(error);
    }

    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    atom_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  // returns the number of cluster at Order=1, which is the number of atoms
  size_t StructureManagerCenters::get_nb_clusters(size_t order) const {
    if (order == 1) {
      return this->n_center_atoms;
    } else {
      throw std::string("ERROR : Order != 1");
    }
  }
  /* ---------------------------------------------------------------------- */

}  // namespace rascal
