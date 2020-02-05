/**
 * file   symmetry_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   17 Dec 2018
 *
 * @brief implementation of symmetry functions for neural nets (G-functions in
 * Behler-Parinello-speak)
 *
 * Copyright © 2018 Till Junge, Markus Stricker COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_REPRESENTATIONS_SYMMETRY_FUNCTIONS_HH_
#define SRC_RASCAL_REPRESENTATIONS_SYMMETRY_FUNCTIONS_HH_

#include "rascal/math/utils.hh"
#include "rascal/structure_managers/property.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/tuple_standardisation.hh"
#include "rascal/utils/units.hh"
#include "rascal/utils/utils.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Eigen/Dense"

namespace rascal {

  using units::UnitStyle;

  enum class SymmetryFunctionType { One, Gaussian, Cosine, Angular1, Angular2 };

  /* ---------------------------------------------------------------------- */
  inline std::string get_name(SymmetryFunctionType fun_type) {
    switch (fun_type) {
    case SymmetryFunctionType::One: {
      return "One";
      break;
    }
    case SymmetryFunctionType::Gaussian: {
      return "Gaussian";
      break;
    }
    case SymmetryFunctionType::Cosine: {
      return "Cosine";
      break;
    }
    case SymmetryFunctionType::Angular1: {
      return "Angular1";
      break;
    }
    case SymmetryFunctionType::Angular2: {
      return "Angular2";
      break;
    }
    default:
      throw std::runtime_error("undefined symmetry function type");
      break;
    }
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType FunType>
  class SymmetryFunction {};

  /* ---------------------------------------------------------------------- */
  std::ostream & operator<<(std::ostream & os,
                            const SymmetryFunctionType & fun_type) {
    os << get_name(fun_type);
    return os;
  }

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType FunType>
  std::ostream & operator<<(std::ostream & os,
                            const SymmetryFunction<FunType> & /*sym_fun*/) {
    os << FunType;
    return os;
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Radial symmetry function as proposed in Behler (2007) PRL 98, 146401; it
   * has two parameters eta and r_s for the form
   * SF(r_ij) = exp(-eta (r_ij - * r_s)^2)
   *
   * where r_ij is the distance between atoms, r_s is the value for shifting and
   * eta scales the exponent.
   */

  template <>
  class SymmetryFunction<SymmetryFunctionType::Gaussian> {
   public:
    static constexpr size_t Order{2};

    using Return_t = std::tuple<double, double>;
    /**
     * usually, derivatives are aligned with the distance vector, in which case
     * a scalar return type is sufficient. (important for triplet-related
     * functions)
     */
    static constexpr bool DerivativeIsCollinear{true};

    // constructor
    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, eta{json_io::check_units(unit_style.distance(-2),
                                                   params.at("eta"))},
          r_s{json_io::check_units(unit_style.distance(), params.at("r_s"))} {}

    double f_sym(const double & r_ij) const {
      auto && delta_r = r_ij - this->r_s;
      return exp(-this->eta * delta_r * delta_r);
    }

    Return_t df_sym(const double & r_ij) const {
      auto && delta_r{r_ij - this->r_s};
      auto && fun_val{exp(-this->eta * delta_r * delta_r)};
      return Return_t(fun_val, -2. * this->eta * delta_r * fun_val);
    }

   protected:
    const json params;
    double eta;
    double r_s;
  };

  constexpr size_t SymmetryFunction<SymmetryFunctionType::Gaussian>::Order;

  /* ---------------------------------------------------------------------- */
  /**
   * Triplet related symmetry function, also called `Angular narrow`. Narrow
   * means that the atom j,k of a triplet also need to be within each others
   * cutoff.
   *
   * SF =
   * 2^(1-zeta) * (1 + lambda * cos_theta)^zeta * exp(r_ij^2 + r_ik^2 + r_jk^2)
   */
  template <>
  class SymmetryFunction<SymmetryFunctionType::Angular1> {
   public:
    static constexpr size_t Order{3};

    // return type to be adjusted?
    using Return_t = std::tuple<double, double>;
    // usage?
    static constexpr bool DerivativeIsCollinear{false};

    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, zeta{json_io::check_units(unit_style.none(),
                                                    params.at("zeta"))},
          lambda{json_io::check_units(unit_style.none(), params.at("lambda"))},
          eta{json_io::check_units(unit_style.distance(-2), params.at("eta"))},
          prefactor{math::pow(2., 1 - zeta)} {}

    double f_sym(const double & cos_theta, const double & r_ij,
                 const double & r_ik, const double & r_jk) const {
      auto && angular_contrib{
          math::pow(1. + this->lambda * cos_theta, this->zeta)};
      auto && exp_contrib{
          exp(-this->eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk))};
      return this->prefactor * angular_contrib * exp_contrib;
    }

    Return_t df_sym(const double & cos_theta, const double & r_ij,
                  const double & r_ik, const double & r_jk) const {
      auto && angular_contrib{
          math::pow(1. + this->lambda * cos_theta, this->zeta)};
      auto && exp_contrib{
          exp(-this->eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk))};
      auto && fun_val{this->prefactor * angular_contrib * exp_contrib};

      // helper for derivative
      // auto && psi_ij{-1/(r_ij * rik) * this->lambda * this->zeta / ()}
      return Return_t(fun_val, 0.);  // placeholder '0' for later
    }

   protected:
    const json params;
    double zeta;
    double lambda;
    double eta;
    double prefactor;  // precomputation at initialization
  };

  // template <>
  // class SymmetryFunction<SymmetryFunctionType::Angular1> {
  //   static constexpr size_t Order{3};
  //   static constexpr size_t NbParams{4};
  //   using ParamShape = Eigen::Matrix<double, NbParams, 1>;
  //   /**
  //    * usually, derivatives are aligned with the distance vector, in which
  //    case
  //    * a scalar return type is sufficient. (important for triplet-related
  //    * functions)
  //    */
  //   static constexpr bool DerivativeIsCollinear{false};

  //   static double eval_function(const Eigen::Map<ParamS> & params,
  //                               const double & r_ij, const double & r_jk,
  //                               const double & r_ik, const double cos_theta)
  //                               {
  //     auto && zeta{params(0)};
  //     auto && eta{params(1)};
  //     auto && lambda{params(2)};
  //     auto && r_s{params(2)};
  //     auto && delta_r = r_ij - r_s;
  //     return exp(-eta * delta_r * delta_r);
  //   }

  //   template <class Derived>
  //   static auto eval_derivative(const Eigen::MatrixBase<ParamShape> & params,
  //                               const double & r_ij,
  //                               const Eigen::MatrixBase<Derived> & n_ij)
  //       -> decltype(auto) {
  //     auto && eta{params(0)};
  //     auto && r_s{params(1)};
  //     auto && delta_r{r_ij - r_s};
  //     return n_ij * (-2. * eta * delta_r * exp(-eta * delta_r * delta_r));
  //   }

  //   static Eigen::Matrix<double, NbParams, 1> read(const json & params,
  //                                                  const UnitStyle & units) {
  //     Eigen::Matrix<double, NbParams, 1> retval{};
  //     retval(0) = json_io::check_units(units.distance(-1, 2),
  //     param.at("eta")); retval(1) = json_io::check_units(units.distance(),
  //     param.at("r_s")); return retval;
  //   }

  //  protected:
  // };
}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_SYMMETRY_FUNCTIONS_HH_
