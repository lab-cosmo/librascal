/**
 * @file   calculator_behler_parrinello_dense.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   12 Sep 2019
 *
 * @brief  Implementation of a dense (i.e. non-sparse) Behler-Parrinello-type
 * descriptor calculater (e.g., as would be normal in a single-species
 * calculation)
 *
 * Copyright © 2019 Till Junge
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "behler_feature.hh"
#include "calculator_base.hh"
#include "symmetry_functions.hh"

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_

namespace rascal {
  // //! forward declaration
  // template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  // class BehlerFeatureBase;

  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  class CalculatorBehlerParrinelloDense : public CalculatorBase {
   public:
    using Parent = CalculatorBase;
    using Hypers_t = typename Parent::Hypers_t;
    using BehlerFeatureBase_t =
        BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>;

    // type of the datastructure used to register the list of valid
    // hyperparameters
    using ReferenceHypers_t = Parent::ReferenceHypers_t;

    //! Default constructor
    CalculatorBehlerParrinelloDense() = delete;

    //! Constructor with input json
    explicit CalculatorBehlerParrinelloDense(
        const Hypers_t & parameters, const UnitStyle & unit_style,
        const std::vector<int> & managed_species_ids,
        const std::map<std::string, int> & species_mapping);

    //! Copy constructor
    CalculatorBehlerParrinelloDense(
        const CalculatorBehlerParrinelloDense & other) = delete;

    //! Move constructor
    CalculatorBehlerParrinelloDense(CalculatorBehlerParrinelloDense && other) =
        default;

    //! Destructor
    virtual ~CalculatorBehlerParrinelloDense() = default;

    //! Copy assignment operator
    CalculatorBehlerParrinelloDense &
    operator=(const CalculatorBehlerParrinelloDense & other) = delete;

    //! Move assignment operator
    CalculatorBehlerParrinelloDense &
    operator=(CalculatorBehlerParrinelloDense && other) = default;

    template <class StructureManager>
    void compute(StructureManager & manager,
                 const std::vector<int> & managed_species_ids) {
      for (auto & behler_feature : this->behler_features) {
        /**
         * todo(jungestricker) BehlerFeatures: what is RepeatedSpecies and
         * permutation, add switch?
         * maybe: throw error if manager does not fit
         */
        if (managed_species_ids.size() != behler_feature->get_order()) {
          // this behler feature handles clusters of a different size and has
          // nothing to do
          break;
        }

        auto expected_repetition{
            get_repeated_species(behler_feature->get_species_ids())};
        switch (expected_repetition) {
        case RepeatedSpecies::Not: {
          this->compute_rep_helper<StructureManager, RepeatedSpecies::Not>(
              manager, *behler_feature, managed_species_ids);
          break;
        }
        case RepeatedSpecies::All: {
          this->compute_rep_helper<StructureManager, RepeatedSpecies::All>(
              manager, *behler_feature, managed_species_ids);
          break;
        }
        case RepeatedSpecies::FirstTwo: {
          this->compute_rep_helper<StructureManager, RepeatedSpecies::FirstTwo>(
              manager, *behler_feature, managed_species_ids);
          break;
        }
        case RepeatedSpecies::SecondTwo: {
          this->compute_rep_helper<StructureManager,
                                   RepeatedSpecies::SecondTwo>(
              manager, *behler_feature, managed_species_ids);
          break;
        }
        case RepeatedSpecies::OuterTwo: {
          this->compute_rep_helper<StructureManager, RepeatedSpecies::OuterTwo>(
              manager, *behler_feature, managed_species_ids);
          break;
        }
        default:
          throw std::runtime_error("Can't handle Unknown RepSpecies");
          break;
        }
      }
    }

    template <class StructureManager, RepeatedSpecies ExpectedRepSpecies>
    void compute_rep_helper(StructureManager & manager,
                            BehlerFeatureBase_t & behler_feature,
                            const std::vector<int> & managed_species_ids) {
      /**
       * todo(tillmarkus): this whole function could be moved to BehlerFeature,
       * because the only thing it does is ask the BehlerFeature about
       * properties and what it should to. Would be nice to have it there.
       * Deferred to later time.
       */
      auto && values{behler_feature.fetch_or_create_value(manager)};
      auto && self_derivatives{
          behler_feature.fetch_or_create_self_derivatives(manager)};
      auto && other_derivatives{
          behler_feature.fetch_or_create_other_derivatives(manager)};

      auto && manager_order{managed_species_ids.size()};
      // the following should be a std::array<int, 4> representing the template
      // parameters of this permutation and should be implemented in
      // permutation.hh
      auto && expected_species_ids{behler_feature.get_species_ids()};
      auto && permutation_label{
          compute_permutation(managed_species_ids, expected_species_ids)};
      switch (manager_order) {
      case PairOrder: {
        switch (permutation_label) {
        case PermutationLabel::p01: {
          using Perm = PermutationSelector<PermutationLabel::p01>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p10: {
          using Perm = PermutationSelector<PermutationLabel::p10>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        default: {
          throw std::runtime_error("unknown Permutation type for pairs");
          break;
        }
        }
        break;
      }
      case TripletOrder: {
        switch (permutation_label) {
        case PermutationLabel::p012: {
          using Perm = PermutationSelector<PermutationLabel::p012>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p120: {
          using Perm = PermutationSelector<PermutationLabel::p120>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p201: {
          using Perm = PermutationSelector<PermutationLabel::p201>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p102: {
          using Perm = PermutationSelector<PermutationLabel::p102>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p021: {
          using Perm = PermutationSelector<PermutationLabel::p021>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        case PermutationLabel::p210: {
          using Perm = PermutationSelector<PermutationLabel::p210>::type;
          behler_feature.template compute<ExpectedRepSpecies, Perm>(
              manager, values, self_derivatives, other_derivatives);
          break;
        }
        default: {
          throw std::runtime_error("unknown Permutation type for triplets");
          break;
        }
        }
        break;
      }

      default: {
        throw std::runtime_error("unknown cluster order");
        break;
      }
      }
    }

   protected:
    /**
     * stores base class refs to the Symmetry functions to be evaluated. The
     * main loop iterates over this container and applies each function to
     * all clusters of an input structure manager, storing the results
     * to a provided input property (which is not necessarily attached to
     * the same structure manager)
     */
    std::vector<std::unique_ptr<BehlerFeatureBase_t>> behler_features{};

    std::string cutoff_function_type;
    std::vector<int> managed_species_ids;
    std::map<std::string, int> species_mapping;
    /**
     * reference of the requiered hypers, these are used to check the
     * parameters for building the behler_features vector at construction
     */
    const ReferenceHypers_t reference_hypers{{"scaling", {}},
                                             {"cutoff_function_type", {}},
                                             {"symmetry_functions", {}}};
  };

  /**
   * Default set of symmetry functions, this is meant to replicate n2p2
   */

  using CalculatorBehlerParrinelloDenseStd =
      CalculatorBehlerParrinelloDense<true, SymmetryFunctionType::Gaussian,
                                      SymmetryFunctionType::AngularNarrow>;
}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_
