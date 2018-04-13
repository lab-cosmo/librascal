/**
 * file   neighbourhood_manager_lammps.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of the neighbourhood manager for lammps
 *        neighbourhood lists
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "neighbourhood_managers/neighbourhood_manager_cell.hh"

#include <numeric>

namespace proteus {

  /* ---------------------------------------------------------------------- */
  /*
  void NeighbourhoodManagerCell::
  reset_impl(const int &inum, const int &tot_num, int *ilist, int *numneigh,
             int **firstneigh, double **x, double **f, int *type,
             double *eatom, double **vatom) {
    this->inum = inum;
    this->tot_num = tot_num;
    this->ilist = ilist;
    this->numneigh = numneigh;
    this->firstneigh = firstneigh;
    this->x = x;
    this->f = f;
    this->type = type;
    this->eatom = eatom;
    this->vatom = vatom;
    this->offsets.reserve(inum);
    this->offsets.resize(1);
    for (int i{0} ; i<this->inum-1 ; ++i) {
      this->offsets.emplace_back(this->offsets[i] + this->numneigh[i]);
    }
    this->nb_pairs = std::accumulate(numneigh, numneigh+this->inum, 0);
  }
  */

  /* ---------------------------------------------------------------------- */
  /*
  size_t NeighbourhoodManagerCell::
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
    default:
      throw std::runtime_error("Can only handle single atoms and pairs");
      break;
    }
  }
  */
}  // proteus
