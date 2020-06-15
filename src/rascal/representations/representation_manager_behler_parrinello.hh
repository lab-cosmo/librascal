/**
 * file   representation_manager_behler_parrinello.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   13 Dec 2018
 *
 * @brief  Behler-Parrinello implementation
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_REPRESENTATIONS_REPRESENTATION_MANAGER_BEHLER_PARRINELLO_HH_
#define SRC_RASCAL_REPRESENTATIONS_REPRESENTATION_MANAGER_BEHLER_PARRINELLO_HH_

#include "cutoff_functions.hh"
#include "representation_manager_base.hh"
#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/property.hh"
#include "structure_managers/species_manager.hh"
#include "utils/tuple_standardisation.hh"

#include <limits>
#include <string>
#include <vector>

namespace rascal {

  template <class StructureManager>
  class BehlerParrinello : public RepresentationManagerBase {
   public:
    static constexpr size_t MaxOrder{StructureManager::traits::MaxOrder};
    static constexpr size_t Dim{StructureManager::traits::Dim};
    static constexpr size_t ForceLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};

    using Parent = RepresentationManagerBase;
    using Force_t = Property<double, 1, ForceLayer, Dim>;
    using StdTypes = TupleStandardisation<MaxOrder>;

    /**
     * structuremanagers are filtered by cutoff for performance reasons
     */
    using StoredStructure_t = AdaptorFilter<StructureManager, MaxOrder>;

    //! Default constructor
    BehlerParrinello() = delete;

    /**
     * Constructor with structure manager and json-formatted hyper parameters
     */
    BehlerParrinello(StructureManager & structure, const json & hypers);

    /**
     * Constructor with structure manager and json-formatted hyper parameters
     */
    BehlerParrinello(StructureManager & structure, const std::string & hypers);

    //! Copy constructor
    BehlerParrinello(const BehlerParrinello & other) = delete;

    //! Move constructor
    BehlerParrinello(BehlerParrinello && other) = default;

    //! Destructor
    virtual ~BehlerParrinello() = default;

    //! Copy assignment operator
    BehlerParrinello & operator=(const BehlerParrinello & other);

    //! Move assignment operator
    BehlerParrinello & operator=(BehlerParrinello && other) = default;

    /**
     * make sure data structures for the compute step are ready. The arguments
     * here are passed up the chain to the structure manager at the root of this
     * stack.
     */
    template <class... Args>
    void update(Args &&... args);

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
    void evaluate_pair_symmetry_function_group(
        std::vector<InputNodeContributionBase> & symmetry_funs);
    constexpr size_t Invalid{std::numeric_limits<size_t>::max()};
    StructureManager & structure;
    SpeciesManager<StructureManager> species;

    //! uninque cutoff function used for all input nodes
    CutoffFuntype cutoff_fun{};
    //! set of all cutoff values for optimisation
    std::set<double> cutoffs{};
    //! number of hidden layers in the network
    /**
     * The next part needs to go out into a kernel-like function which only
     * deals with the neural network
     */
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

    std::map<StdTypes, std::vector<InputNodeContributionBase>>
        symmetry_functions;

    const std::string symmetry_function_key{"symmetry function values"};
    const std::string symmetry_derivative_key{"symmetry function derivatives"};
  };

}  // namespace rascal

#include "representation_manager_behler_parrinello_impl.hh"

#endif  // SRC_RASCAL_REPRESENTATIONS_REPRESENTATION_MANAGER_BEHLER_PARRINELLO_HH_
