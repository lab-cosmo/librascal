/**
 * file   structure_manager_json.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   30 May 2018
 *
 * @brief implementation structure manager for reading atomic structure
 *        from file
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/structure_manager_json.hh"

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
  void StructureManagerJson::update() {

    /**
     * Check if a structure has already been read, if not, throw an exception to
     * let the user know, what is wrong.
     */
    try {
      auto pos_size = atoms_object.position.size();
      if (pos_size == 0) {
        throw std::runtime_error("No atomic structure defined. "
                                 "Read structure first!");
      }
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    for (auto vec : atoms_object.cell) {
      for (auto coord : vec) {
        this->cell_data.push_back(coord);
      }
    }

    for (auto pos : atoms_object.position) {
      for (auto coord : pos) {
        this->pos_data.push_back(coord);
      }
    }

    /**
     * Before going further, we check if the number of positions and atom types
     * match, otherwise the data set is incomplete.
     */
    assert(atoms_object.position.size() == atoms_object.type.size());
    /**
     * Set the protected member variable number of atoms, depending on the given
     * number of positions
     */
    this->natoms = atoms_object.position.size();
    /**
     * The following two commands build a list of increasing indices. It is
     * assumed that the atoms do not have a unique index, when they are read
     * from file. Therefore a list of increasing integer identifiers is
     * built. Numbers are assigned to the positions in the order in which they
     * appear in the file.
     */
    this->ilist.resize(this->natoms);
    std::iota(ilist.begin(), ilist.end(), 0);

    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};

    atom_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Helper function to get the number of <code>clusters</code> (atoms, pairs,
   * triplets, quadruplets) with a specific <code>Order</code>. The
   * <code>MaxOrder</code> depends on your implementation or
   * processing. Increasing the order is done with an adaptor.
   */
  size_t StructureManagerJson::get_nb_clusters(size_t cluster_size) const {
    switch (cluster_size) {
    case 1: {
      return this->natoms;
      break;
    }
    default:
      throw std::runtime_error("Can only handle atoms and pairs,"
                               " use adaptor to increase MaxOrder.");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Function for reading data from a JSON file in the ASE format. See the
   * definition of <code>AtomicStructure</code> and adapt the fields, which
   * should be read to your case. One peculiarity should be mentioned: The type
   * <code>std::fstream</code> does not throw an exception, if the file is not
   * read. That is the reason for the try/catch block -- to make sure, the file
   * is opened.
   */
  void StructureManagerJson::
  read_structure_from_json(const std::string filename) {

    json j;

    try {
      std::ifstream f(filename);
      if (!f.is_open()) throw std::ios::failure("Error opening JSON file!");
      f >> j;
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    /**
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
    this->atoms_object = j.begin().value();
  }
} // rascal
