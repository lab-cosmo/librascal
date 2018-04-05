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
 * @section LICENSE
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

#include "neighbourhood_managers/neighbourhood_manager_lammps.hh"

#include <numeric>

namespace proteus {

  /* ---------------------------------------------------------------------- */
  NeighbourhoodManagerLammps::
  NeighbourhoodManagerLammps(const int & inum, const int & tot_num,
                             int * ilist, int * numneigh, int ** firstneigh,
                             double ** x, double ** f, int * type,
                             double * eatom, double ** vatom)
    : inum{inum}, tot_num{tot_num}, ilist{ilist}, numneigh{numneigh},
      firstneigh{firstneigh}, x{x}, f{f}, type{type}, eatom{eatom},
      vatom{vatom}, nb_pairs{std::accumulate(numneigh, numneigh+inum, 0)}
  { }

}  // proteus
