/**
 * file   neighbourhood_manager_chain.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   30 May 2018
 *
 * @brief Implementation of the neighbourhood manager for polyalanine
 *        chain from json file
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

#include "neighbourhood_managers/neighbourhood_manager_chain.hh"

#include <numeric>
#include <fstream>
#include <iostream>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  void NeighbourhoodManagerChain::update() {
    // Ensure contiguous data structures
    for (const auto vec : neigh_in.cell) {
      for (const auto coord : vec) {
	this->cell_data.push_back(coord);
      }
    }

    for (const auto pos : neigh_in.position) {
      for (const auto coord : pos) {
	this->pos_data.push_back(coord);
      }
    }

    this->natoms = neigh_in.position.size();
  }

  /* ---------------------------------------------------------------------- */
  size_t NeighbourhoodManagerChain::get_nb_clusters(int cluster_size)  {
    switch (cluster_size) {
    case 1: {
      return natoms;
      break;
    }
    default:
      throw std::runtime_error("Can only handle single atoms; "
                               " use adaptor to increase MaxLevel.");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  void NeighbourhoodManagerChain::
  read_structure_from_json(const std::string filename) {

    json j;

    try {
      std::ifstream f(filename);
      if (!f.is_open()) throw std::ios::failure("Error opening JSON file!");
      f >> j;
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
    }

    // ASE json format is nested - here, first entry is actual molecule
    this->neigh_in = j.begin().value();
  }

  /* ---------------------------------------------------------------------- */


} // rascal
