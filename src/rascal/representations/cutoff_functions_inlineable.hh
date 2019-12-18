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

  class CutoffFunctionBase {
   public:
    //! Constructor
    CutoffFunctionBase(const InlCutoffFunctionType & cut_fun_type)
        : cut_fun_type{cut_fun_type} {}
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
                        const bool & compute_derivatives) const;

    /**
     * The identifier string must provide a unique name for a property to
     * store precomputed cutoff function values. For this, the naming scheme
     * must guarantee that two different parameter sets (e.g., cut-off radii)
     * generate diffent names, and should guarantee that the same function
     * with the same  parameters generates the same name.
     */
    virtual const std::string & get_identifier() const = 0;

   protected:
    InlCutoffFunctionType cut_fun_type;
    //! Main worker (raison d'être)
    template <InlCutoffFunctionType CutFunType, class StructureManager>
    inline void compute_helper(StructureManager & manager,
                               const bool & compute_derivatives) const;
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
        : Parent{InlCutoffFunctionType::Cosine}, hypers{hypers} {
      // make sure the hypers contain a parameter r_cut
      cutoff = json_io::check_units(unit_style.distance(), hypers.at("r_cut"));
      identifier = this->make_identifier(unit_style);
    }

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
    Hypers_t hypers;
    //! cutoff radii
    double cutoff{};
    std::string identifier{};
  };

  template <class StructureManager>
  void CutoffFunctionBase::compute(StructureManager & manager,
                                   const bool & compute_derivatives) const {
    switch (this->cut_fun_type) {
    case InlCutoffFunctionType::Cosine: {
      this->template compute_helper<InlCutoffFunctionType::Cosine>(
          manager, compute_derivatives);
      break;
    }
    default:
      throw std::runtime_error("Unknown cutoff function type");
      break;
    }
  }

  template <InlCutoffFunctionType CutFunType, class StructureManager>
  inline void
  CutoffFunctionBase::compute_helper(StructureManager & manager,
                                     const bool & compute_derivatives) const {
    auto & typed_this{static_cast<const CutoffFunction<CutFunType> &>(*this)};

    // check whether the property already exists (assuming everyone computes
    // either values or both values and derivatives
    const std::string value_identifier{this->get_identifier() + "_value"};
    const std::string derivative_identifier{this->get_identifier() +
                                            "_derivative"};
    auto & identifier{compute_derivatives ? derivative_identifier
                                          : value_identifier};

    bool property_at_wrong_level{
        manager.is_property_in_stack(identifier) and
        not manager.is_property_in_current_level(identifier)};
    if (property_at_wrong_level) {
      // complain and die
      throw std::runtime_error(
          "cannot handle the situation where the property isn't registered "
          "at the same stack level as this cutoff function");
    }

    using Property_t = Property<double, PairOrder, StructureManager>;
    constexpr bool Validate{false}, AllowCreation{true};

    auto & property{*manager.template get_property<Property_t>(
        identifier, Validate, AllowCreation)};

    if (property.is_updated()) {
      return;
    }

    property.resize();

    // compute cutoff functions
    if (not compute_derivatives) {
      for (auto && atom : manager) {
        for (auto && pair : atom.pairs()) {
          property[pair] = typed_this.f_c(manager.get_distance(pair));
        }
      }
    } else {
      auto & value_property{*manager.template get_property<Property_t>(
          value_identifier, Validate, AllowCreation)};
      value_property.resize();

      for (auto && atom : manager) {
        for (auto && pair : atom.pairs()) {
          std::tie(value_property[pair], property[pair]) =
            std::tuple_cat(typed_this.df_c(manager.get_distance(pair)));
        }
      }
    }
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_INLINEABLE_HH_
