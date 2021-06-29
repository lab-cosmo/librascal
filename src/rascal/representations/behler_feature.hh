/**
 * file   behler_feature.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief utilities for the definition of input nodes (combinations of symmetry
 * functions with cutoff functions
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_
#define SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_

#include "rascal/representations/cutoff_functions_inlineable.hh"
#include "rascal/representations/symmetry_functions.hh"
#include "rascal/structure_managers/property.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/permutation.hh"

namespace rascal {

  /**
   * A BehlerFeature is a single G function with a single set of parameters
   * @tparam CompatibilityMode Behler's and Singraber's implementations evaluate
   * triplets where the j and k atom have the same species only once, despite
   * using full neighbour lists. This is not a problem mathematically, but
   * inconsistent with the mathematical description in the papers.
   * CompatibilityMode = true replicates the inconsistent behaviour for
   * convenience in comparisons.
   */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  class BehlerFeatureBase {
   public:
    constexpr static int MaxBehlerOrder{3};
    using StdSpecies = TupleStandardisation<int, MaxBehlerOrder>;

    using Hypers_t = json;

    template <class StructureManager>
    using GVals_t = Property<double, AtomOrder,
                             StructureManager>;
    template <class StructureManager>
    using dGSelfVals_t =
        Property<double, AtomOrder,
                 StructureManager, ThreeD>;

    template <class StructureManager>
    using dGOtherVals_t =
        Property<double, PairOrder,
                 StructureManager, ThreeD, 2>;

    constexpr static bool CompatibilityMode{CompatibilityMode_};
    //! Default constructor
    BehlerFeatureBase() = delete;

    //! Constructor with symmetry function type
    BehlerFeatureBase(const SymmetryFunctionType & sym_fun_type,
                      std::shared_ptr<CutoffFunctionBase> cut_fun,
                      const size_t & order, const Hypers_t & raw_params)
        : sym_fun_type{sym_fun_type}, order{order},
          raw_params{raw_params}, cut_fun{cut_fun} {}

    /**
     * factory function for unique ptr from json
     */
    inline static std::unique_ptr<BehlerFeatureBase>
    make_unique(std::shared_ptr<CutoffFunctionBase> cut_fun,
                const UnitStyle & unit_style, const Hypers_t & parameters,
                const std::map<std::string, int> & species_numbers);

    //! Copy constructor
    BehlerFeatureBase(const BehlerFeatureBase & other);

    //! Move constructor
    BehlerFeatureBase(BehlerFeatureBase && other) noexcept;

    //! Destructor
    virtual ~BehlerFeatureBase() = default;

    //! Copy assignment operator
    BehlerFeatureBase & operator=(const BehlerFeatureBase & other);

    //! Move assignment operator
    BehlerFeatureBase & operator=(BehlerFeatureBase && other) noexcept;

    //! needs to be called after reading the input file and prior to the first
    //! evaluation. Attaches all necessary properties for precalculated values
    //! to the manager
    virtual void init(const UnitStyle & units) = 0;

    //! Main worker (raison d'être) computes input node values
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    inline void compute(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output_values) const;

    //! Main worker (raison d'être) computes input node values and derivatives
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    inline void
    compute(StructureManager & manager,
            std::shared_ptr<PropertyBase> output_values,
            std::shared_ptr<PropertyBase> output_self_derivatives,
            std::shared_ptr<PropertyBase> output_other_derivatives) const;

    //! insert a parameter (sub-)json
    void add_params(const json & params) {
      const auto & type{json_io::get<std::string>(params, "type")};
      if (type != get_name(this->sym_fun_type)) {
        std::stringstream error{};
        error << "Parameter set for function type '" << type
              << "' assigned to function of type '"
              << get_name(this->sym_fun_type) << "'.";
        throw std::runtime_error(error.str());
      }
      this->raw_params.push_back(params);
    }

    /**
     * returns a string fully qualifying the feature
     */
    std::string get_identifier() const {
      std::stringstream name_stream{};
      name_stream << "Behlerfeature_cutfun_" << this->cut_fun->get_identifier()
                  << "_symfun_" << this->get_sym_fun_identifier();
      return name_stream.str();
    }

    template <class StructureManager>
    std::shared_ptr<PropertyBase>
    fetch_or_create_value(StructureManager & manager) {
      constexpr bool ValidateProperty{true};
      constexpr bool AllowCreation{true};
      return manager.template get_property<GVals_t<StructureManager>>(
          this->get_identifier() + "_value", ValidateProperty, AllowCreation);
    }

    template <class StructureManager>
    std::shared_ptr<PropertyBase>
    fetch_or_create_self_derivatives(StructureManager & manager) {
      constexpr bool ValidateProperty{true};
      constexpr bool AllowCreation{true};
      return manager.template get_property<dGSelfVals_t<StructureManager>>(
          this->get_identifier() + "_self_derivatives", ValidateProperty,
          AllowCreation);
    }

    template <class StructureManager>
    std::shared_ptr<PropertyBase>
    fetch_or_create_other_derivatives(StructureManager & manager) {
      constexpr bool ValidateProperty{true};
      constexpr bool AllowCreation{true};
      return manager.template get_property<dGOtherVals_t<StructureManager>>(
          this->get_identifier() + "_other_derivatives", ValidateProperty,
          AllowCreation);
    }

    virtual const std::vector<int> & get_species_ids() const = 0;
    const size_t & get_order() const { return this->order; }

   protected:
    virtual std::string get_sym_fun_identifier() const = 0;
    template <SymmetryFunctionType... FunTypes_>
    struct SymFunctionsVTable;

    const SymmetryFunctionType sym_fun_type;
    const size_t order;
    std::vector<json> raw_params{};
    RepeatedSpecies species_repetition{RepeatedSpecies::Unknown};
    // StdSpecies species_combo{};
    bool is_initialised{false};
    std::shared_ptr<CutoffFunctionBase> cut_fun;
  };

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  constexpr int
      BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>::MaxBehlerOrder;
  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  constexpr bool
      BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>::CompatibilityMode;

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  class BehlerPairFeature final
      : public BehlerFeatureBase<CompatibilityMode_, SymFunTypes...> {
   public:
    using SymmetryFunction_t = SymmetryFunction<MySymFunType>;
    // Pair features don't need a compatibility mode, see desctription above
    using Parent = BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>;
    using Hypers_t = typename Parent::Hypers_t;

    static constexpr size_t Order{SymmetryFunction_t::Order};
    static_assert(Order == PairOrder,
                  "Should only be instantiated for pair symmetry functions");

    //! Default constructor
    BehlerPairFeature(std::shared_ptr<CutoffFunctionBase> cut_fun,
                      const UnitStyle & unit_style, const Hypers_t & raw_params,
                      const std::map<std::string, int> & species_numbers)
        : Parent(MySymFunType, cut_fun, SymmetryFunction_t::Order, raw_params),
          sym_fun{
              unit_style,
              json_io::get(raw_params, "params" + canary(raw_params, "params")),
              species_numbers} {
      if (json_io::get<std::string>(raw_params, "type") !=
          get_name(MySymFunType)) {
        std::stringstream err{};
        err << "params for symmetry function of type '"
            << json_io::get<std::string>(raw_params, "type")
            << "' provided to initialise a symmetry function of type '"
            << get_name(MySymFunType);
        throw std::runtime_error(err.str());
      }

      if (cut_fun->get_cutoff() !=
          json_io::check_units(
              unit_style.distance(),
              json_io::get(json_io::get(raw_params, "cutoff_function"),
                           "r_cut"))) {
        std::stringstream err{};
        err << "Mismatch: the provided cutoff function has a cuttoff radius of "
            << cut_fun->get_cutoff()
            << " but the parameters prescribe a cutoff radius of "
            << json_io::get<double>(
                   json_io::get(json_io::get(raw_params, "cutoff_function"),
                                "r_cut"),
                   "value");
        throw std::runtime_error(err.str());
      }
    }

    //! Copy constructor
    BehlerPairFeature(const BehlerPairFeature & other) = delete;

    //! Move constructor
    BehlerPairFeature(BehlerPairFeature && other) = default;

    //! Destructor
    ~BehlerPairFeature() = default;

    //! Copy assignment operator
    BehlerPairFeature & operator=(const BehlerPairFeature & other) = delete;

    //! Move assignment operator
    BehlerPairFeature & operator=(BehlerPairFeature && other) = default;

    void init(const UnitStyle & units) final;

    const std::vector<int> & get_species_ids() const final {
      return this->sym_fun.get_species_ids();
    }

    // template <RepeatedSpecies RepSpecies, class StructureManager>
    // void compute(StructureManager & manager,
    //              std::shared_ptr<PropertyBase> output) const;

    // template <RepeatedSpecies RepSpecies, class StructureManager>
    // void compute(StructureManager & manager,
    //              std::shared_ptr<PropertyBase> output,
    //              std::shared_ptr<PropertyBase> output_derivatives) const;

    size_t get_index() const;  // to implement

    /* ---------------------------------------------------------------------- */
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    void compute_helper(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output) const;

    /* ---------------------------------------------------------------------- */
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    void compute_helper(
        StructureManager & manager, std::shared_ptr<PropertyBase> output,
        std::shared_ptr<PropertyBase> output_self_derivatives,
        std::shared_ptr<PropertyBase> output_other_derivatives) const;

   protected:
    std::string get_sym_fun_identifier() const final {
      return this->sym_fun.get_identifier();
    }
    SymmetryFunction<MySymFunType> sym_fun;
  };

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  constexpr size_t BehlerPairFeature<CompatibilityMode_, MySymFunType,
                                     SymFunTypes...>::Order;

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode, SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  class BehlerTripletFeature final
      : public BehlerFeatureBase<CompatibilityMode, SymFunTypes...> {
   public:
    using SymmetryFunction_t = SymmetryFunction<MySymFunType>;
    using Parent = BehlerFeatureBase<CompatibilityMode, SymFunTypes...>;
    using Hypers_t = typename Parent::Hypers_t;

    static constexpr size_t Order{SymmetryFunction_t::Order};
    static_assert(Order == TripletOrder,
                  "Should only be instantiated for triplet symmetry functions");
    //! Default constructor
    BehlerTripletFeature(std::shared_ptr<CutoffFunctionBase> cut_fun,
                         const UnitStyle & unit_style,
                         const Hypers_t & raw_params,
                         const std::map<std::string, int> & species_numbers)
        : Parent(MySymFunType, cut_fun, SymmetryFunction_t::Order, raw_params),
          sym_fun{unit_style, json_io::get(raw_params, "params"),
                  species_numbers} {
      if (json_io::get<std::string>(raw_params, "type") !=
          get_name(MySymFunType)) {
        std::stringstream err{};
        err << "params for symmetry function of type '"
            << json_io::get<std::string>(raw_params, "type")
            << "' provided to initialise a symmetry function of type '"
            << get_name(MySymFunType);
        throw std::runtime_error(err.str());
      }

      if (cut_fun->get_cutoff() !=
          json_io::check_units(
              unit_style.distance(),
              json_io::get(json_io::get(raw_params, "cutoff_function"),
                           "r_cut"))) {
        std::stringstream err{};
        err << "Mismatch: the provided cutoff function has a cuttoff radius of "
            << cut_fun->get_cutoff()
            << " but the parameters prescribe a cutoff radius of "
            << json_io::get_quantity(raw_params, "r_cut");
        throw std::runtime_error(err.str());
      }
    }

    //! Copy constructor
    BehlerTripletFeature(const BehlerTripletFeature & other) = delete;

    //! Move constructor
    BehlerTripletFeature(BehlerTripletFeature && other) = default;

    //! Destructor
    ~BehlerTripletFeature() = default;

    //! Copy assignment operator
    BehlerTripletFeature &
    operator=(const BehlerTripletFeature & other) = delete;

    //! Move assignment operator
    BehlerTripletFeature & operator=(BehlerTripletFeature && other) = default;

    void init(const UnitStyle & units) final;

    const std::vector<int> & get_species_ids() const final {
      return this->sym_fun.get_species_ids();
    }

    // template <RepeatedSpecies RepSpecies, class StructureManager>
    // void compute(StructureManager & manager,
    //              std::shared_ptr<PropertyBase> output) const;

    // template <RepeatedSpecies RepSpecies, class StructureManager>
    // void compute(StructureManager & manager,
    //              std::shared_ptr<PropertyBase> output,
    //              std::shared_ptr<PropertyBase> output_derivatives) const;

    size_t get_index() const;  // to implement

    /* ---------------------------------------------------------------------- */
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    void compute_helper(StructureManager & manager,
                        std::shared_ptr<PropertyBase> output) const;

    /* ---------------------------------------------------------------------- */
    template <RepeatedSpecies RepSpecies, typename Permutation,
              class StructureManager>
    void compute_helper(
        StructureManager & manager, std::shared_ptr<PropertyBase> output,
        std::shared_ptr<PropertyBase> output_self_derivatives,
        std::shared_ptr<PropertyBase> output_other_derivatives) const;

   protected:
    std::string get_sym_fun_identifier() const final {
      return this->sym_fun.get_identifier();
    }
    SymmetryFunction<MySymFunType> sym_fun;
  };

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  constexpr size_t BehlerTripletFeature<CompatibilityMode_, MySymFunType,
                                        SymFunTypes...>::Order;

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  std::unique_ptr<BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>>
  BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>::make_unique(
      std::shared_ptr<CutoffFunctionBase> cut_fun, const UnitStyle & unit_style,
      const Hypers_t & parameters,
      const std::map<std::string, int> & species_numbers) {
    auto && sym_fun_label{json_io::get<std::string>(parameters, "type")};
    if (sym_fun_label == "Gaussian") {
      return std::make_unique<BehlerPairFeature<
          CompatibilityMode_, SymmetryFunctionType::Gaussian, SymFunTypes...>>(
          cut_fun, unit_style, parameters, species_numbers);
    } else if (sym_fun_label == "AngularNarrow") {
      return std::make_unique<BehlerTripletFeature<
          CompatibilityMode_, SymmetryFunctionType::AngularNarrow,
          SymFunTypes...>>(cut_fun, unit_style, parameters, species_numbers);
    } else if (sym_fun_label == "AngularWide") {
      return std::make_unique<BehlerTripletFeature<
          CompatibilityMode_, SymmetryFunctionType::AngularWide,
          SymFunTypes...>>(cut_fun, unit_style, parameters, species_numbers);
    } else {
      throw std::runtime_error("Unknown symmetry function type '" +
                               sym_fun_label + "'");
    }
  }

}  // namespace rascal
#include "behler_feature_impl.hh"

#endif  // SRC_RASCAL_REPRESENTATIONS_BEHLER_FEATURE_HH_
