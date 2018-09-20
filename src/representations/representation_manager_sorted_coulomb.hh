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
#include <map>
#include <vector>

namespace rascal {
  // //! forward declaration for traits
  // class RepresentationManagerSortedCoulomb;

  // template <>
  // struct RepresentationManager_traits<RepresentationManagerSortedCoulomb>
  // {
  //   constexpr static int Dim{3};
  //   constexpr static size_t MaxOrder{1};
  //   constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::yes};
  //   //using LayerByDimension = std::integer_sequence<size_t, 0, 0>;
  // };
  template<class StructureManager>
  class RepresentationManagerSortedCoulomb//: public RepresentationManagerBase
  {
  public:
    // using traits = RepresentationManager_traits<RepresentationManagerSortedCoulomb>;
    using hypers_t = RepresentationManagerBase::hypers_t;
    using Manager_t = StructureManager;
    //! Default constructor 
    // RepresentationManagerSortedCoulomb()
    //   :structure_manager{},central_decay{},interaction_cutoff{},interaction_decay{}
    //   {}
    RepresentationManagerSortedCoulomb(StructureManager &sm, 
      double central_decay , double interaction_cutoff, 
      double interaction_decay)
      :structure_manager{sm},central_decay{central_decay},
      interaction_cutoff{interaction_cutoff},interaction_decay{interaction_decay}
      {}
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
    
    // Pure Virtual Function to set hyperparameters of the representation
    // TODO think of a generic input type for the hypers
    // virtual void set_hyperparameters(const hypers_t & );

    void build(StructureManager &sm, 
      double central_decay , double interaction_cutoff, 
      double interaction_decay){
        this->structure_manager = sm;
        this->central_decay = central_decay;
        this->interaction_cutoff = interaction_cutoff;
        this->interaction_decay = interaction_decay;
      };

    Manager_t& structure_manager;
    //hypers_t hyperparmeters;
    double central_decay;
    double interaction_cutoff;
    double interaction_decay;

  protected:
  private:
  };

}

#endif /* BASIS_REPRESENTATION_MANAGER_SORTED_COULOMB_H */

