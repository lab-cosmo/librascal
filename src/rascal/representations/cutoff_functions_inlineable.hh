/**
 * @file   cutoff_functions_inlineable.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   09 Dec 2019
 *
 * @brief  bundle inlineable cutoff functions
 *
 * Copyright © 2019 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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
 *
 */

#ifndef SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_INLINEABLE_HH_
#define SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_INLINEABLE_HH_

#include "rascal/representations/calculator_base.hh"
#include "rascal/structure_managers/property.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/units.hh"

#include <iostream>

namespace rascal {

  /**
   * List of implemented cutoff function
   */
  enum class InlCutoffFunctionType { Cosine, End_ };

  /**
   * forward declaration for compute dispatch
   */
  template <InlCutoffFunctionType Type>
  class CutoffFunction;

  enum class Evaluation { Value, Derivative };
  class CutoffFunctionBase {
   public:
    template <typename StructureManager>
    using Property_t = Property<double, PairOrder, StructureManager>;

    //! Constructor
    explicit CutoffFunctionBase(const InlCutoffFunctionType & cut_fun_type,
                                const double & cutoff)
        : cut_fun_type{cut_fun_type}, cutoff{cutoff} {}
    //! Destructor
    virtual ~CutoffFunctionBase() = default;
    //! Copy constructor
    CutoffFunctionBase(const CutoffFunctionBase & other) = delete;
    //! Move constructor
    CutoffFunctionBase(CutoffFunctionBase && other) = default;
    //! Copy assignment operator
    CutoffFunctionBase & operator=(const CutoffFunctionBase & other) = delete;
    //! Move assignment operator
    CutoffFunctionBase & operator=(CutoffFunctionBase && other) = default;

    using Hypers_t = CalculatorBase::Hypers_t;

    //! Main worker (raison d'être)
    template <class StructureManager>
    inline void compute(StructureManager & manager,
                        const Evaluation & evaluation) const;

    /**
     * The identifier string must provide a unique name for a property to
     * store precomputed cutoff function values. For this, the naming scheme
     * must guarantee that two different parameter sets (e.g., cut-off radii)
     * generate diffent names, and should guarantee that the same function
     * with the same  parameters generates the same name.
     */
    virtual const std::string & get_identifier() const = 0;

    /**
     * returns the property with cutoff functions value, guaranteed to be
     * fresh
     */
    template <typename StructureManager>
    inline Property_t<StructureManager> &
    get_value(StructureManager & manager) const;

    /**
     * returns the property with cutoff functions values and derivatives,
     * guaranteed to be fresh
     */
    template <typename StructureManager>
    std::tuple<
        Property_t<StructureManager> &,
        Property_t<StructureManager> &> inline get_derivative(StructureManager &
                                                                  manager)
        const;

    const double & get_cutoff() const { return this->cutoff; }

   protected:
    template <typename StructureManager>
    Property_t<StructureManager> &
    get_property(StructureManager & manager,
                 const Evaluation & evaluation) const;
    InlCutoffFunctionType cut_fun_type;
    //! Main worker (raison d'être)
    template <InlCutoffFunctionType CutFunType, class StructureManager>
    inline void compute_helper(StructureManager & manager,
                               const Evaluation & evaluation) const;
    //! cutoff radii
    const double cutoff;
  };

  /**
   * Cosine cutoff function as in Behler, can only be used with strict
   * managers
   */
  template <>
  class CutoffFunction<InlCutoffFunctionType::Cosine>
      : public CutoffFunctionBase {
   public:
    using Parent = CutoffFunctionBase;
    using Hypers_t = typename CutoffFunctionBase::Hypers_t;
    explicit CutoffFunction(const units::UnitStyle & unit_style,
                            const Hypers_t & hypers)
        : Parent{InlCutoffFunctionType::Cosine,
                 json_io::check_units(unit_style.distance(),
                                      hypers.at("r_cut"))},
          hypers{hypers}, identifier{this->make_identifier(unit_style)} {}

    inline double f_c(const double & distance) const {
      assert(distance <= this->cutoff);
      return .5 * (std::cos(math::PI * distance / this->cutoff) + 1.);
    }

    inline std::array<double, 2> df_c(const double & distance) const {
      assert(distance <= this->cutoff);
      auto && scaled_dist{math::PI * distance / this->cutoff};
      return {.5 * (std::cos(scaled_dist) + 1.),
              -.5 * scaled_dist * std::sin(scaled_dist)};
    }

    const std::string & get_identifier() const { return this->identifier; }

   protected:
    std::string make_identifier(const units::UnitStyle & unit_style) {
      std::stringstream id{};
      id.precision(14);
      id << "Cosine_" << this->cutoff << '_' << unit_style.distance();
      return id.str();
    }
    //! keep the hypers
    const Hypers_t hypers;
    const std::string identifier;
  };

  template <class StructureManager>
  void CutoffFunctionBase::compute(StructureManager & manager,
                                   const Evaluation & evaluation) const {
    switch (this->cut_fun_type) {
    case InlCutoffFunctionType::Cosine: {
      this->template compute_helper<InlCutoffFunctionType::Cosine>(manager,
                                                                   evaluation);
      break;
    }
    default:
      throw std::runtime_error("Unknown cutoff function type");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  template <typename StructureManager>
  auto CutoffFunctionBase::get_value(StructureManager & manager) const
      -> Property_t<StructureManager> & {
    constexpr Evaluation EvalKind{Evaluation::Value};
    this->compute(manager, EvalKind);
    return this->get_property(manager, EvalKind);
  };

  /* ---------------------------------------------------------------------- */
  template <typename StructureManager>
  auto CutoffFunctionBase::get_derivative(StructureManager & manager) const
      -> std::tuple<Property_t<StructureManager> &,
                    Property_t<StructureManager> &> {
    this->compute(manager, Evaluation::Derivative);
    return std::make_tuple(this->get_property(manager, Evaluation::Value),
                           this->get_property(manager, Evaluation::Derivative));
  };

  /* ---------------------------------------------------------------------- */
  template <typename StructureManager>
  auto CutoffFunctionBase::get_property(StructureManager & manager,
                                        const Evaluation & evaluation) const
      -> Property_t<StructureManager> & {
    // check whether the property already exists (assuming everyone computes
    // either values or both values and derivatives)
    const std::string value_identifier{this->get_identifier() + "_value"};
    const std::string derivative_identifier{this->get_identifier() +
                                            "_derivative"};
    auto & identifier{evaluation == Evaluation::Derivative
                          ? derivative_identifier
                          : value_identifier};

    bool property_at_wrong_layer{
        manager.is_property_in_stack(identifier) and
        not manager.is_property_in_current_level(identifier)};
    if (property_at_wrong_layer) {
      // complain and die
      throw std::runtime_error(
          "cannot handle the situation where the property isn't registered "
          "at the same stack layer as this cutoff function");
    }

    constexpr bool Validate{false}, AllowCreation{true};

    auto & property{
        *manager.template get_property<Property_t<StructureManager>>(
            identifier, Validate, AllowCreation)};

    return property;
  }

  /* ---------------------------------------------------------------------- */
  template <InlCutoffFunctionType CutFunType, class StructureManager>
  inline void
  CutoffFunctionBase::compute_helper(StructureManager & manager,
                                     const Evaluation & evaluation) const {
    auto & typed_this{static_cast<const CutoffFunction<CutFunType> &>(*this)};
    auto & property{this->get_property(manager, evaluation)};
    if (property.is_updated()) {
      return;
    }

    property.resize();

    switch (evaluation) {
    case Evaluation::Value: {
      for (auto && atom : manager) {
        for (auto && pair : atom.pairs()) {
          property[pair] = typed_this.f_c(manager.get_distance(pair));
        }
      }
      property.set_updated_status(true);
      break;
    }
    case Evaluation::Derivative: {
      auto & value_property{this->get_value(manager)};
      value_property.resize();

      for (auto && atom : manager) {
        for (auto && pair : atom.pairs()) {
          std::tie(value_property[pair], property[pair]) =
              std::tuple_cat(typed_this.df_c(manager.get_distance(pair)));
        }
      }
      property.set_updated_status(true);
      value_property.set_updated_status(true);
      break;
    }
    default:
      throw std::runtime_error("Unknown evaluation type");
      break;
    }

  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_INLINEABLE_HH_
