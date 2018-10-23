/**
 * file   representation_manager_sorted_coulomb.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  base class for representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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


#include "representations/representation_manager_sorted_coulomb.hh"

namespace rascal {

  template<class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::set_hyperparameters(
          const RepresentationManagerSortedCoulomb<Mngr>::hypers_t & hyper){

    this->central_decay = hyper["central_decay"];
    this->interaction_cutoff = hyper["interaction_cutoff"];
    this->interaction_decay = hyper["interaction_decay"];
    this->size = hyper["size"];
  };

  template<class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::build(){
    
    // upper diag of the coulomb mat
    Eigen::MatrixXd lin_coulomb{};
    lin_coulomb.resize(this->size*(this->size+1)/2);

    // upper diag of the coulomb mat
    std::vector<double> distances_to_sort{};
    //! initialise the coulomb_matrices storage
    this->coulomb_matrices.resize_to_zero();

    for (auto center: this->structure_manager){

      for (auto neigh1: center){
        int ii{neigh1.get_index()};
        auto distance{this->structure_manager.get_distance(neigh1)};
        distances_to_sort.push_back(distance);

        for (auto neigh2: center){
          int jj{neigh2.get_index()};
          if (ii >= jj) continue;

          auto dij{(neigh1.get_position()-neigh2.get_position()).norm};
          
        }
      }
    }
  }

}