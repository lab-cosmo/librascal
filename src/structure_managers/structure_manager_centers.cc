/**
 * file   structure_manager_centers.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   06 August 2018
 *
 * @brief Manager with atoms and centers
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "structure_managers/structure_manager_centers.hh"

#include <numeric>
#include <fstream>
#include <iostream>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  /**
   * After reading the <code>atoms_object</code> from the file, the cell vectors
   * as well as the atomic positions are put into contiguous a
   * <code>std::vector</code> data structure to ensure fast access via the
   * <code>Eigen::Map</code>.  <code>std::vector</code> provide iterator access,
   * which is used here with the <code>auto</code> keyword. Using this, it is
   * unnecessary to know the exact size of the positions/cell array. No
   * distinction between 2 or 3 dimensions. We just put all numbers in a vector
   * and access them with the map. Using the vector type automatically ensures
   * contiguity
   */

  void StructureManagerCenters::
  update(const Eigen::Ref<const Eigen::MatrixXd, 0,
         Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> positions,
         const Eigen::Ref<const VecXi> atom_types,
         const Eigen::Ref<const Eigen::MatrixXd> cell,
         const Eigen::Ref<const PBC_t> pbc) {

    Eigen::Index Natom{positions.cols()};
    this->natoms = Natom;
    this->atoms_object.set_structure(positions, atom_types, cell, pbc);

    StructureManagerCenters::build();
    
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    atom_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  void StructureManagerCenters::build() {

    this->atoms_object.atoms_type.resize(this->natoms);
    //set the references to the particles positions
    for (size_t id{0}; id < this->natoms; ++id){
      this->atoms_index[0].push_back(id);
      this->offsets.push_back(id);
    }
    Cell_t lat = this->atoms_object.cell;
    this->lattice.set_cell(lat);
  }

  /* ---------------------------------------------------------------------- */
  size_t StructureManagerCenters::get_nb_clusters(size_t /*order*/)
    const {
    return this->natoms;
  }

} // rascal
