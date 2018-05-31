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
#include "external/json.hpp"

#include <numeric>





namespace rascal {

  /* ---------------------------------------------------------------------- */
  void NeighbourhoodManagerChain::
  update() {
    // create neighbourlist/triplet etc.
  }


  /* ---------------------------------------------------------------------- */
  size_t NeighbourhoodManagerChain::
  get_nb_clusters(int cluster_size)  {
    switch (cluster_size) {
    case 1: {
      return inum;
      break;
    }
    case 2: {
      return nb_pairs;
      break;
    }
    case 3: {
      return nb_triplets;
      break;
    }
    default:
      throw std::runtime_error("Can only handle single atoms,"
                               " pairs and triplets");
      break;
    }
  }

}  // rascal
