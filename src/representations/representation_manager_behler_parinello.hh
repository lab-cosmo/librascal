/**
 * file   representation_manager_behler_parinello.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   13 Dec 2018
 *
 * @brief  Behler-Parinello implementation
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with librascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "representation_manager_base.hh"
#include "cutoff_functions.hh"
#include "structure_managers/species_manager.hh"
#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/property.hh"

#include <string>
#include <vector>
#include <limits>

#ifndef REPRESENTATION_MANAGER_BEHLER_PARINELLO_H
#define REPRESENTATION_MANAGER_BEHLER_PARINELLO_H

namespace rascal {

  template <class StructureManager>
  class BehlerParinello : public RepresentationManagerBase {
   public:
    static constexpr size_t MaxOrder{StructureManager::traits::MaxOrder};
    static constexpr size_t Dim{StructureManager::traits::Dim};
    static constexpr size_t ForceLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};

    using Parent = RepresentationManagerBase;
    using Force_t = Property<double, 1, ForceLayer, Dim>;

    /**
     * structuremanagers are filtered by cutoff for performance reasons
     */
    using StoredStructure_t = AdaptorFilter<StructureManager, MaxOrder>;

    //! Default constructor
    BehlerParinello() = delete;

    /**
     * Constructor with structure manager and json-formatted hyper parameters
     */
    BehlerParinello(StructureManager & structure, const json & hypers);

    /**
     * Constructor with structure manager and json-formatted hyper parameters
     */
    BehlerParinello(StructureManager & structure, const std::string & hypers);

    //! Copy constructor
    BehlerParinello(const BehlerParinello & other) = delete;

    //! Move constructor
    BehlerParinello(BehlerParinello && other) = default;

    //! Destructor
    virtual ~BehlerParinello() = default;

    //! Copy assignment operator
    BehlerParinello & operator=(const BehlerParinello & other);

    //! Move assignment operator
    BehlerParinello & operator=(BehlerParinello && other) = default;

    /**
     * Evaluate all features
     */
    void compute() final;

    /**
     * Evaluate all features
     */
    void evaluate_forces();

    //! Pure Virtual Function to set hyperparameters of the representation
    void set_hyperparameters(const Hypers_t &) final;

    //! get the raw data of the representation
    std::vector<Precision_t> & get_representation_raw_data() final {
      throw std::runtime_error("does not apply");
    }

    //! get the size of a feature vector
    size_t get_feature_size() final {
      throw std::runtime_error("does not apply");
    }

    //! get the number of centers for the representation
    size_t get_center_size() final { return this->structure.size(); }

   protected:
    constexpr size_t Invalid{std::numeric_limits<size_t>::max()};
    StructureManager & structure;
    SpeciesManager<StructureManager> species;

    // neural net topology
    //! uninque cutoff function used for all input nodes
    CutoffFuntype cutoff_fun{};
    //! set of all cutoff values for optimisation
    std::set<double> cutoffs{};
    //! number of hidden layers in the network
    size_t nb_hidden_layers{Invalid};
    //! number of species represented by the network (some species might be
    //! missing at runtime)
    size_t nb_species{Invalid};
    //! number of nodes per hidden layer, counting in direction of output coming
    //! from input
    std::vector<size_t> nb_per_hidden{Invalid};
    //! per species: number of symmetry functions defined with that same species
    //! defined as primary
    std::map<int, size_t> nb_sym_per_species{Invalid};

    std::map<double, SymmetryFunction> sym
    
    
    std::map<double,
             typename SpeciesManager<StructureManager>::FilterContainer_t>
        manager_by_cutoff;
  };

}  // namespace rascal

#include "representation_manager_behler_parinello_impl.hh"

#endif /* REPRESENTATION_MANAGER_BEHLER_PARINELLO_H */
