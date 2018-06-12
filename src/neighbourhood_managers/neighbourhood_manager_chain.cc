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
    this->make_neighbourlist();
  }

  /* ---------------------------------------------------------------------- */
  size_t NeighbourhoodManagerChain::get_nb_clusters(int cluster_size)  {
    switch (cluster_size) {
    case 1: {
      return this->natoms;
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

    // ASE json format is nested - here, first entry is actual data structure
    this->neigh_in = j.begin().value();
  }

  /* ---------------------------------------------------------------------- */
  inline double NeighbourhoodManagerChain::get_box_length(int d) {
    Cell_ref Cell = this->get_cell();
    return Cell.col(d).norm();
  }
  /* ---------------------------------------------------------------------- */
  inline int get_linear_index(std::vector<int> nidx, std::vector<int> nmax) {
    auto dim = nidx.size();
    switch (dim) {
    case 1: {
      return nidx[0];
      break;
    }
    case 2: {
      return nidx[1]*nmax[0] + nidx[0];
      break;
    }
    case 3: {
      return nidx[2] * nmax[0] * nmax[1] +
	nidx[1] * nmax[0] + nidx[0];
      break;
    }
    default:
      throw std::runtime_error("Can only give index for 1,2,3 dimensions");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  void NeighbourhoodManagerChain::make_neighbourlist() {
    // internal variables for linked list/ linked cell
    std::vector<int> nmax(3);
    std::vector<double> rc(3);

    for(auto dim{0}; dim < traits::Dim; ++dim) {
      nmax[dim] = static_cast<int>(std::floor(this->get_box_length(dim)
					    / this->cut_off));
      rc[dim] = static_cast<double>(this->get_box_length(dim) / nmax[dim]);
    };

    int nboxes{1};
    for (auto n : nmax){nboxes *= n;}
    nboxes = std::max(nboxes, 1);

    std::cout << "nboxes " << nboxes << std::endl;
    std::cout << "nmax "
	      << nmax[0] << " "
	      << nmax[1] << " "
	      << nmax[2] << " "
	      << std::endl;
    std::vector<int> ll(this->natoms);
    std::vector<int> lc(nboxes) ;
    for (auto & i : ll){
      i = -1;
    }
    for (int & i : lc) {
      i = -1;
    }

    Positions_ref atom_pos = this->get_positions();
    Eigen::Matrix<double, 1, traits::Dim> offset{};
    for(auto dim{0}; dim < traits::Dim; ++dim) {
      offset(dim) = atom_pos.row(dim).minCoeff();
      std::cout << "box length " << this->get_box_length(dim) << std::endl;
    }
    std::cout<< "offset " << offset << std::endl;

    // Make cell lists
    std::vector<int> nidx(traits::Dim);
    for (auto i{0}; i < atom_pos.cols(); ++i) {
      auto * p{atom_pos.col(i).data()};
      Vector_ref pos{p};
      std::cout << "p " << pos << std::endl;

      for(auto dim{0}; dim < traits::Dim; ++dim) {
	// std::cout << "rc[dim] " << " "
	// 	  << pos(dim)-offset(dim)
	// 	  << " " <<  rc[dim] << std::endl;
	nidx[dim] = static_cast<int>(std::floor( (pos(dim) - offset(dim))
						 / rc[dim] ));
	// std::cout << "nidx[dim] 1 " << nidx[dim] << std::endl;
	nidx[dim] = std::min(nidx[dim], nmax[dim] - 1);
	// std::cout << "nidx[dim] 2 " << nidx[dim] << std::endl;
	nidx[dim] = std::max(nidx[dim], 0);
	// std::cout << "nidx[dim] 3 " << nidx[dim] << std::endl;
      }

      auto linear_index = get_linear_index(nidx, nmax);

      std::cout<< "Linear index: " << linear_index << std::endl;
      std::cout<< "dim-index: "
	       << nidx[0] << " "
	       << nidx[1] << " "
	       << nidx[2] << " "
	       << std::endl;

      ll[i] = lc[linear_index];
      lc[linear_index] = i;

    }

    // print neighbour list
    std::cout << ">>>> nboxes " << nboxes << std::endl;
    for (auto i{0}; i<nboxes; ++i) {
      auto n = lc[i];
      std::cout << "linear index " << i << std::endl;
      while (n != -1) {
	std::cout << " " << n << std::endl;
	n = ll[n];
      }
    }

    // Make actual neighbourlist arrays





  }


} // rascal
