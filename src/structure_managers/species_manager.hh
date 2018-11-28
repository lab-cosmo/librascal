/**
 * file   species_manager.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   04 May 2018
 *
 * @brief iterable proxy to a neigbourhood, filtered by species. You
 * can iterate over it by atom species, and get a subset of the
 * neighbourhood of only pairs, triplets ... for which the first atom
 * is of a given type and descend recursively
 *
 * @section LICENSE
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */



#ifndef SPECIES_MANAGER_H
#define SPECIES_MANAGER_H

#include "structure_managers/structure_manager.hh"
#include "representations/basis_function_manager.hh"

namespace rascal {

  template <class NeighManager, int MaxLayer, int Layer = 0>
  class SpeciesManager {
   public:
    using Species_t = int;
    using Dummy_t = char;
    using SubManagerType =
      std::conditional_t<(Layer != MaxLayer),
                         SpeciesManager<NeighManager, MaxLayer, Layer+1>,
                         Dummy_t>;
    //! Default constructor
    SpeciesManager() = delete;

    //! Construct from an existing StructureManager
    explicit SpeciesManager(NeighManager & manager);

    //! Copy constructor
    SpeciesManager(const SpeciesManager &other) = delete;

    //! Move constructor
    SpeciesManager(SpeciesManager &&other) = default;

    //! Destructor
    virtual ~SpeciesManager()  = default;

    //! Copy assignment operator
    SpeciesManager& operator=(const SpeciesManager &other) = delete;

    //! Move assignment operator
    SpeciesManager& operator=(SpeciesManager &&other) = default;

    //! get the symmetry functions and the corresponding StructureManager
    std::map<Species_t, BasisFunManager> & get_symmetry_functions();

    //! get the next depth layer
    template <bool NotAtMaxLayer = (Layer != MaxLayer)>
    std::map<Species_t,
             std::enable_if_t<NotAtMaxLayer,
                              SubManagerType>>
    & get_next_order();

   protected:
    std::array<Species_t, Layer> fixed_species;
    std::map<Species_t, BasisFunManager>;
    std::map<Species_t, SubManagerType>;
  };

}  // rascal

#endif /* SPECIES_MANAGER_H */
