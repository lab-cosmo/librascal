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
    // Invoke neighbourlist, because without it, you can not really to anything
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
  inline std::vector<int>
  NeighbourhoodManagerChain::
  get_box_index(Vector_ref& position,
		std::vector<double>& rc,
		Eigen::Matrix<double, 1, traits::Dim>offset,
		std::vector<int> nmax) {
    std::vector<int> nidx(traits::Dim);
    for(auto dim{0}; dim < traits::Dim; ++dim) {
      nidx[dim] = static_cast<int>(std::floor( (position(dim) - offset(dim))
					       / rc[dim] ));
      nidx[dim] = std::min(nidx[dim], nmax[dim] - 1);
      nidx[dim] = std::max(nidx[dim], 0);
    }
    return nidx;
  }

  /* ---------------------------------------------------------------------- */
  inline void NeighbourhoodManagerChain::
  collect_neighbour_info_of_atom(const int i,
				  const std::vector<int> boxidx,
				  const std::vector<int> nmax) {

    auto jcell_index = get_linear_index(boxidx, nmax);
    auto ihead = this->lc[jcell_index]; //ll[jcell_index];
    // Check if any atom is in given cell
    if(ihead != -1) {
      if(ihead != i) this->firstneigh[i].push_back(ihead);
      // same cell neighbours?
      auto itail = this->ll[ihead];
      while (itail != -1) {
	if(itail != i) {
	  this->firstneigh[i].push_back(itail);
	}
	itail = this->ll[itail];
      }
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

    this->ll.resize(this->natoms);
    this->lc.resize(nboxes) ;
    for (auto & i : this->ll){
      i = -1;
    }
    for (int & i : this->lc) {
      i = -1;
    }

    Positions_ref atom_pos = this->get_positions();
    Eigen::Matrix<double, 1, traits::Dim> offset{};
    for(auto dim{0}; dim < traits::Dim; ++dim) {
      offset(dim) = std::min(0., atom_pos.row(dim).minCoeff());
    }

    // Make cell lists
    std::vector<int> nidx(traits::Dim);
    for (auto i{0}; i < atom_pos.cols(); ++i) {
      auto * p{atom_pos.col(i).data()};
      Vector_ref pos{p};
      nidx = get_box_index(pos, rc, offset, nmax);
      auto linear_index = get_linear_index(nidx, nmax);
      this->ll[i] = this->lc[linear_index];
      this->lc[linear_index] = i;

    }

    // print cell list
    std::cout << ">>>> nboxes " << nboxes << std::endl;
    for (auto i{0}; i<nboxes; ++i) {
      auto n = this->lc[i];
      // std::cout << "linear index " << i << std::endl;
      while (n != -1) {
	std::cout << "box " << i<< " atom " << n << std::endl;
	n = this->ll[n];
      }
    }
    // Make verlet table? of neighbours (vectorized)
    // Full neighbour list
    this->firstneigh.resize(this->natoms);
    std::cout << ">>>>>>>>>>>>>>>>start neighbour list" << std::endl;
    for (auto i{0}; i < atom_pos.cols(); ++i) {
      auto * p{atom_pos.col(i).data()};
      Vector_ref pos{p};
      // auto ineigh{0};
      std::cout << "---pos \n" << pos << std::endl;

      // own cell
      nidx = get_box_index(pos, rc, offset, nmax);

      auto nneighcells {std::pow(3, traits::Dim) - 1};
      std::cout << "Number of neighbour cells " << nneighcells << std::endl;

      auto nidxtmp = nidx;
      switch(traits::Dim) {
      case 1: {
	std::vector<int> nd{-1, 0, 1};
	for (auto dx :nd) {
	  nidxtmp[0] = nidx[0] + dx ;
	  if(nidxtmp[0] < 0) continue;
	  if(nidxtmp[1] > nmax[1] -1) continue;
	  collect_neighbour_info_of_atom(i, nidxtmp, nmax);
	}
	break;
      }
      case 2: {
	std::vector<int> nd{-1, 0, 1};
	for (auto dx : nd) {
	  for (auto dy : nd) {
	    nidxtmp[0] = nidx[0] + dx ;
	    nidxtmp[1] = nidx[1] + dy ;
	    if(nidxtmp[0] < 0) continue;
	    if(nidxtmp[0] > nmax[0] -1) continue;
	    if(nidxtmp[1] < 0) continue;
	    if(nidxtmp[1] > nmax[1] -1) continue;
	    collect_neighbour_info_of_atom(i, nidxtmp, nmax);
	  }
	}
	break;
      }
      case 3: {
	// Full neighbourlist
	std::vector<int> nd{-1, 0, 1};
	// Half neighbour list?
	// std::vector<int> nd{0, 1};
	for (auto dx : nd) {
	  for (auto dy : nd) {
	    for (auto dz : nd) {
	      nidxtmp[0] = nidx[0] + dx ;
	      nidxtmp[1] = nidx[1] + dy ;
	      nidxtmp[2] = nidx[2] + dz ;
	      if(nidxtmp[0] < 0) continue;
	      if(nidxtmp[0] > nmax[0] -1) continue;
	      if(nidxtmp[1] < 0) continue;
	      if(nidxtmp[1] > nmax[1] -1) continue;
	      if(nidxtmp[2] < 0) continue;
	      if(nidxtmp[2] > nmax[2] -1) continue;
	      collect_neighbour_info_of_atom(i, nidxtmp, nmax);
	    }
	  }
	}
	break;
      }
      default:
	throw std::runtime_error("Neighbourlist only max 3D.");
	break;
      }


      // loop through neighbours
    } // atom neighbours

    // for (auto a : firstneigh) {
    //   std::cout << "===" << std::endl;
    //   for (auto b : a) {
    // 	// b++;
    // 	std::cout << "neighbours "<< b << std::endl;
    //   }
    //   // std::cout << std::endl;
    // }
    // Half-neighbourlist
    this->numneigh.resize(this->natoms);
    this->nb_pairs = 0;

    for (size_t i=0; i < this->natoms; ++i) {
      auto neigh = firstneigh[i];
      auto i_number_of_neighbours{0};
      for(size_t j : neigh) {
	if (j > i) {
	  i_number_of_neighbours += 1;
	  this->nb_pairs += 1;
	  halfneigh.push_back(j);
	}
	numneigh[i] = i_number_of_neighbours;
      }
    }

    std::cout << "Number of pairs: " << this->nb_pairs << std::endl;
    auto ioff{0};
    for (size_t n=0; n < this->natoms; ++n){
      auto a = numneigh[n];
      std::cout << "Next neigh atom---- " << n
		<< " ynumber of neighbours " << a << std::endl;
      for (auto i=ioff; i < ioff+a; i++) {
	std::cout << "Neighbour "
		  << halfneigh[i] << std::endl;
      }
      ioff += a;
    }







  }


} // rascal
