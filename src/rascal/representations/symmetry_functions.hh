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

  enum class SymmetryFunctionType { Gaussian, AngularNarrow, AngularWide };

  /* ---------------------------------------------------------------------- */
  inline std::string get_name(SymmetryFunctionType fun_type) {
    switch (fun_type) {
    case SymmetryFunctionType::Gaussian: {
      return "Gaussian";
      break;
    }
    case SymmetryFunctionType::AngularNarrow: {
      return "AngularNarrow";
      break;
    }
    case SymmetryFunctionType::AngularWide: {
      return "AngularWide";
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
   * means that the atoms j,k of a triplet also need to be within each others
   * cutoff.
   *
   * SF =
   * 2^(1-zeta) * (1 + lambda * cos_theta)^zeta * exp(r_ij^2 + r_ik^2 + r_jk^2)
   */
  template <>
  class SymmetryFunction<SymmetryFunctionType::AngularNarrow> {
   public:
    static constexpr size_t Order{3};

    // return type for each function value and 3 derivative values
    using Return_t = std::tuple<double, double, double, double>;
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

      double dval_i{0};
      double dval_j{0};
      double dval_k{0};
      return Return_t(fun_val, dval_i, dval_j, dval_k);
    }

   protected:
    const json params;
    double zeta;
    double lambda;
    double eta;
    double prefactor;  // prefactor evaluation at construction
  };

  constexpr size_t SymmetryFunction<SymmetryFunctionType::AngularNarrow>::Order;

  /* ---------------------------------------------------------------------- */
  /**
   * Triplet related symmetry function, also called `Angular wide`. Wide means
   * that the atoms j,k of a triplet do not need to be within each others
   * cutoff.
   *
   * SF =
   * 2^(1-zeta) * (1 + lambda * cos_theta)^zeta * exp(r_ij^2 + r_ik^2)
   */
  template <>
  class SymmetryFunction<SymmetryFunctionType::AngularWide> {
   public:
    static constexpr int DistsPerTriplet{3};
    static constexpr size_t Order{3};

    // return type for each function value and 3 derivative values
    using Return_t = std::tuple<double, std::array<double, DistsPerTriplet>>;
    // usage?
    static constexpr bool DerivativeIsCollinear{false};

    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, zeta{json_io::check_units(unit_style.none(),
                                                    params.at("zeta"))},
          lambda{json_io::check_units(unit_style.none(), params.at("lambda"))},
          eta{json_io::check_units(unit_style.distance(-2), params.at("eta"))},
          prefactor{math::pow(2., 1 - zeta)} {}

    double f_sym(const double & cos_theta, const double & r_ij,
                 const double & r_ik, const double & /*r_jk*/) const {
      auto && angular_contrib{
          math::pow(1. + this->lambda * cos_theta, this->zeta)};
      auto && exp_contrib{exp(-this->eta * (r_ij * r_ij + r_ik * r_ik))};
      return this->prefactor * angular_contrib * exp_contrib;
    }

    Return_t df_sym(const double & cos_theta, const double & r_ij,
                    const double & r_ik, const double & r_jk) const {
      auto && lam_cos_theta_1 {1. + this->lambda * cos_theta};
      auto && angular_contrib{
          math::pow(lam_cos_theta_1, this->zeta)};
      auto && r_ij2{r_ij * r_ij};
      auto && r_ik2{r_ik * r_ik};
      auto && exp_contrib{exp(-this->eta * (r_ij2 + r_ik2))};
      auto && fun_val{this->prefactor * angular_contrib * exp_contrib};

      auto && d_dr_ij{-2 * this->eta * r_ij * fun_val +
                      this->zeta * (this->lambda / r_ik - cos_theta / r_ij) *
                          fun_val / lam_cos_theta_1};
      auto && d_dr_ik{-2 * this->eta * r_ik * fun_val +
                      this->zeta * (this->lambda / r_ij - cos_theta / r_ik) *
                          fun_val / lam_cos_theta_1};
      auto && d_dr_jk{-this->lambda * r_jk * this->zeta * fun_val /
                      (r_ij * r_ik * lam_cos_theta_1)};
      return Return_t(fun_val, {d_dr_ij, d_dr_ik, d_dr_jk});
    }

   protected:
    const json params;
    double zeta;
    double lambda;
    double eta;
    double prefactor;  // prefactor (2^(1-ζ))
  };

  constexpr size_t SymmetryFunction<SymmetryFunctionType::AngularWide>::Order;

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_SYMMETRY_FUNCTIONS_HH_
