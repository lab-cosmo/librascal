/**
 * file   representation_manager_sorted_coulomb.hh
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

#ifndef BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H
#define BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H

#include "representations/representation_manager_base.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include <vector>

namespace rascal {


  template<class StructureManager>
  class RepresentationManagerSortedCoulomb: public RepresentationManagerBase
  {
  public:
    // TODO make a traits mechanism
    using hypers_t = RepresentationManagerBase::hypers_t;
    using Property_t = Property<double, 1, 1, Eigen::Dynamic, 1>;
    using Manager_t = StructureManager;

    //! Default constructor 
    RepresentationManagerSortedCoulomb(Manager_t &sm, 
      double central_decay , double interaction_cutoff, 
      double interaction_decay, size_t size)
      :structure_manager{sm},central_decay{central_decay},
      interaction_cutoff{interaction_cutoff},
      interaction_decay{interaction_decay},size{size},
      coulomb_matrices{sm}
      {
        // this->compute();
      }

    RepresentationManagerSortedCoulomb(Manager_t &sm,const hypers_t& hyper)
      :structure_manager{sm},central_decay{},
      interaction_cutoff{},
      interaction_decay{},coulomb_matrices{sm}
      {
        this->set_hyperparameters(hyper);
      }

    //! Copy constructor
    RepresentationManagerSortedCoulomb(
      const RepresentationManagerSortedCoulomb &other) = delete;

    //! Move constructor
    RepresentationManagerSortedCoulomb(
      RepresentationManagerSortedCoulomb &&other) = default;

    //! Destructor
    virtual ~RepresentationManagerSortedCoulomb()  = default;

    //! Copy assignment operator
    RepresentationManagerSortedCoulomb& operator=(
      const RepresentationManagerSortedCoulomb &other) = delete;

    //! Move assignment operator
    RepresentationManagerSortedCoulomb& operator=(
      RepresentationManagerSortedCoulomb && other) = default;
    

    // get representation
    void compute();


    // TODO think of a generic input type for the hypers
    void set_hyperparameters(const hypers_t & );


    Manager_t& structure_manager;
    //hypers_t hyperparmeters;
    double central_decay;
    double interaction_cutoff;
    double interaction_decay;
    // first dimension of the largest coulomb mat
    size_t size;

    Property_t coulomb_matrices;

  protected:
  private:
  };


  template<class Mngr>
  void RepresentationManagerSortedCoulomb<Mngr>::compute(){
    
    // upper diag of the coulomb mat
    Eigen::MatrixXd lin_coulomb(this->size,this->size);
    // lin_coulomb.resize(this->size*(this->size+1)/2);
    // lin_coulomb.resize(this->size*this->size);

    // upper diag of the coulomb mat
    std::vector<double> distances_to_sort{};
    //! initialise the coulomb_matrices storage
    this->coulomb_matrices.resize_to_zero();

    for (auto center: this->structure_manager){

      for (auto neigh_i: center){
        size_t ii{neigh_i.get_index()};
        auto dik{this->structure_manager.get_distance(neigh_i)};
        distances_to_sort.push_back(dik);
        auto Zi{neigh_i.get_atom_type()};
        lin_coulomb(ii*this->size) = 0.5*std::pow(Zi,2.4);
        for (auto neigh_j: center){
          size_t jj{neigh_j.get_index()};
          // work only on the lower diagonal
          if (ii > jj) continue;
          auto Zj{neigh_i.get_atom_type()};
          auto dij{(neigh_i.get_position()-neigh_j.get_position()).norm()};
          lin_coulomb(ii*this->size+jj) = Zi*Zj/dij;
        }
      }
    }
  }

}

#endif /* BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H */

