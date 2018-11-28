/**
 * file   structure_manager_centers.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   06 August 2018
 *
 * @brief Manager with atoms and centers
 *
 * Copyright © 2018  Felix Musil, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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
  /*
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
         const Eigen::Ref<const Eigen::VectorXi> atom_types,
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
  // overloading the update function to be able to update from a file
  void StructureManagerCenters::update(const std::string filename) {
    this->read_structure_from_json(filename);

    this->natoms = this->atoms_object.positions.size() / traits::Dim;
    StructureManagerCenters::build();

    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    atom_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  // function for setting the internal data structures
  void StructureManagerCenters::build() {
    // set the references to the particles positions
    for (size_t id{0}; id < this->natoms; ++id) {
      this->atoms_index[0].push_back(id);
      this->offsets.push_back(id);
    }

    Cell_t lat = this->atoms_object.cell;
    this->lattice.set_cell(lat);
  }

  /* ---------------------------------------------------------------------- */
  // returns the number of cluster at Order=1, which is the number of atoms
  size_t StructureManagerCenters::get_nb_clusters(size_t order) const {
    if (order == 1) {
      return this->natoms;
    } else{
      throw std::string("ERREUR : Order != 1");
    }
  }

  /* ---------------------------------------------------------------------- */
  /*
   * One peculiarity should be mentioned: The type <code>std::fstream</code>
   * does not throw an exception, if the file is not read. That is the reason
   * for the try/catch block -- to make sure, the file is opened.
   */
  void StructureManagerCenters::
  read_structure_from_json(const std::string filename) {
    // atoms object do hold read data from a file
    json_io::AtomicJsonData json_atoms_object{};

    json j;

    try {
      std::ifstream f(filename);
      if (!f.is_open()) throw std::ios::failure("Error opening JSON file!");
      f >> j;
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    /*
     * ASE json format is nested - here, first entry is actual data
     * structure. If in any case you should have multiple
     * <code>atoms_objects</code> in your file, which you want to read, the
     * following line has to be adapted. Nesting on the first level is
     * structure1, structure 2, etc. These could be a time series of a
     * simulation, but also just different structures you want to read in from
     * different simulation runs.  Each structure should contain the necessary
     * fields for the <code>AtomicStructure</code> object defined in the header
     * belonging to this file. Here, just the first one is read.
     */
    json_atoms_object = j.begin().value();
    this->atoms_object.set_structure(json_atoms_object);
  }

  /* ---------------------------------------------------------------------- */

} // rascal
