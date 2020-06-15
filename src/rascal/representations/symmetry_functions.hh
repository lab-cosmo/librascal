/**
 * file   symmetry_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   17 Dec 2018
 *
 * @brief implementation of symmetry functions for neural nets (G-functions in
 * Behler-Parrinello-speak)
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
  class SymmetryFunctionBase {
   public:
    using PairReturn_t = std::tuple<double, double>;
    using TripletReturn_t =
        std::tuple<double, std::array<double, nb_distances(TripletOrder)>>;

    double f_sym(const double & /*dummy*/) const {
      throw std::runtime_error("wrong order");
    }
    template <class Derived0, class Derived1, class Derived2>
    double f_sym(const Eigen::MatrixBase<Derived0> & /*cos_theta*/,
                 const Eigen::MatrixBase<Derived1> & /*distances*/,
                 const Eigen::MatrixBase<Derived2> & /*cutoff_vals*/,
                 const std::array<size_t, 3> & /*ordering*/ = {0, 1, 2}) const {
      throw std::runtime_error("wrong order");
    }

    PairReturn_t df_sym(const double & /*dummy*/) const {
      throw std::runtime_error("wrong order");
    }
    template <class Derived0, class Derived1, class Derived2, class Derived3>
    TripletReturn_t
    df_sym(const Eigen::MatrixBase<Derived0> & cos_theta,
           const Eigen::MatrixBase<Derived1> & distances,
           const Eigen::MatrixBase<Derived2> & cutoff_vals,
           const Eigen::MatrixBase<Derived3> & cutoff_derivatives,
           const std::array<size_t, 3> & ordering = {0, 1, 2}) const {
      throw std::runtime_error("wrong order");
    }
  };
  /* ---------------------------------------------------------------------- */
  template <SymmetryFunctionType FunType>
  class SymmetryFunction : public SymmetryFunctionBase {};

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
   * SF(r_ij) = exp(-η (r_ij - * r_s)²)
   *
   * where r_ij is the distance between atoms, r_s is the value for shifting and
   * eta scales the exponent.
   */

  std::string canary(const json & params, const std::string & name) {
    try {
      std::cout << params.at(name) << std::endl;
    } catch (json::exception & error) {
      std::stringstream errmsg{};
      errmsg << "Can't access field '" << name << "' in json '" << params
             << "'";
      throw std::runtime_error(errmsg.str());
    }
    return "";
  }

  template <>
  class SymmetryFunction<SymmetryFunctionType::Gaussian>
      : public SymmetryFunctionBase {
   public:
    using Parent = SymmetryFunctionBase;
    using Parent::df_sym;
    using Parent::f_sym;

    static constexpr size_t Order{2};

    using Return_t = Parent::PairReturn_t;
    /**
     * usually, derivatives are aligned with the distance vector, in which case
     * a scalar return type is sufficient. (important for triplet-related
     * functions)
     */
    static constexpr bool DerivativeIsCollinear{true};

    // constructor
    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, eta{json_io::check_units(
                              unit_style.distance(-2),
                              params.at("eta" + canary(params, "eta")))},
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
   * 2^(1-ζ) * (1 + λ * cos(ϑ))^ζ * exp(r_ij² + r_ik² + r_jk²)
   */
  template <>
  class SymmetryFunction<SymmetryFunctionType::AngularNarrow>
      : public SymmetryFunctionBase {
   public:
    static constexpr size_t Order{3};

    using Parent = SymmetryFunctionBase;
    using Parent::df_sym;
    using Parent::f_sym;

    // return type for each function value and 3 derivative values
    using Return_t = Parent::TripletReturn_t;
    // usage?
    static constexpr bool DerivativeIsCollinear{false};

    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, zeta{json_io::check_units(unit_style.none(),
                                                    params.at("zeta"))},
          lambda{json_io::check_units(unit_style.none(), params.at("lambda"))},
          eta{json_io::check_units(unit_style.distance(-2), params.at("eta"))},
          prefactor{math::pow(2., 1 - zeta)} {}

    template <class Derived0, class Derived1, class Derived2>
    double f_sym(const Eigen::MatrixBase<Derived0> & cos_theta,
                 const Eigen::MatrixBase<Derived1> & distances,
                 const Eigen::MatrixBase<Derived2> & cutoff_vals,
                 const std::array<size_t, 3> & ordering = {0, 1, 2}) const {
      static_assert(Derived0::SizeAtCompileTime == nb_distances(Order),
                    "Distance cos_theta has wrong size");
      static_assert(Derived1::SizeAtCompileTime == nb_distances(Order),
                    "Distance array has wrong size");
      static_assert(Derived2::SizeAtCompileTime == nb_distances(Order),
                    "Cutoff function values array has wrong size");
      auto && angular_contrib{
          math::pow(1. + this->lambda * cos_theta[ordering[0]], this->zeta)};
      auto && exp_contrib{exp(-this->eta * distances.squaredNorm())};
      return this->prefactor * angular_contrib * exp_contrib *
             cutoff_vals.prod();
    }

    template <class Derived0, class Derived1, class Derived2, class Derived3>
    Return_t df_sym(const Eigen::MatrixBase<Derived0> & cos_thetas,
                    const Eigen::MatrixBase<Derived1> & distances,
                    const Eigen::MatrixBase<Derived2> & /*cutoff_vals*/,
                    const Eigen::MatrixBase<Derived3> & /*cutoff_derivatives*/,
                    const std::array<size_t, 3> & ordering = {0, 1, 2}) const {
      auto && cos_theta{cos_thetas(0)};
      auto && lam_cos_theta_1{1. + this->lambda * cos_theta};
      auto && angular_contrib{math::pow(lam_cos_theta_1, this->zeta)};
      auto && exp_contrib{exp(-this->eta * distances.squaredNorm())};
      auto && fun_val{this->prefactor * angular_contrib * exp_contrib};

      // todo(jungestricker) order of distances, should be fine according to
      // Singraber (2019) J Chem Theory Comput 15, 1827; the ordering of
      // triplet distances in the argument here is {d_ij, d_jk, d_ki}, i.e. a
      // circular ordering. But the equations are written as {r_1, r_2, r_3}
      // which correspond to our r_1 -> r_ij, r_2 -> r_ki, r_3 -> r_jk

      auto && r_ij{distances(ordering[0])};
      auto && r_jk{distances(ordering[1])};
      auto && r_ki{distances(ordering[2])};

      double d_dr_ij{
          -2 * this->eta * r_ij * fun_val +
          this->zeta * (this->lambda / r_ki - this->lambda / r_ij * cos_theta) *
              fun_val / lam_cos_theta_1};
      double d_dr_ki{
          -2 * this->eta * r_ki * fun_val +
          this->zeta * (this->lambda / r_ij - this->lambda / r_ki * cos_theta) *
              fun_val / lam_cos_theta_1};
      double d_dr_jk{-this->lambda * r_jk * this->zeta * fun_val /
                         (r_ij * r_ki * lam_cos_theta_1) -
                     2 * this->eta * r_jk * fun_val};

      return Return_t(fun_val, {d_dr_ij, d_dr_jk, d_dr_ki});
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
   * 2^(1-ζ) * (1 + λ * cos(ϑ))^ζ * exp(r_ij² + r_ik²)
   */
  template <>
  class SymmetryFunction<SymmetryFunctionType::AngularWide>
      : public SymmetryFunctionBase {
   public:
    using Parent = SymmetryFunctionBase;
    using Parent::df_sym;
    using Parent::f_sym;

    static constexpr size_t Order{3};

    // return type for each function value and 3 derivative values
    using Return_t = Parent::TripletReturn_t;
    // usage?
    static constexpr bool DerivativeIsCollinear{false};

    SymmetryFunction(const UnitStyle & unit_style, const json & params)
        : params{params}, zeta{json_io::check_units(unit_style.none(),
                                                    params.at("zeta"))},
          lambda{json_io::check_units(unit_style.none(), params.at("lambda"))},
          eta{json_io::check_units(unit_style.distance(-2), params.at("eta"))},
          prefactor{math::pow(2., 1 - zeta)} {}

    template <class Derived0, class Derived1, class Derived2>
    double f_sym(const Eigen::MatrixBase<Derived0> & cos_theta,
                 const Eigen::MatrixBase<Derived1> & distances,
                 const Eigen::MatrixBase<Derived2> & /*cutoff_vals*/,
                 const std::array<size_t, 3> & ordering = {0, 1, 2}) const {
      auto && angular_contrib{
          math::pow(1. + this->lambda * cos_theta[ordering[0]], this->zeta)};
      auto && r_ij{distances[ordering[0]]};
      auto && r_ki{distances[ordering[2]]};
      auto && exp_contrib{exp(-this->eta * (r_ij * r_ij + r_ki * r_ki))};
      return this->prefactor * angular_contrib * exp_contrib;
    }

    template <class Derived0, class Derived1, class Derived2, class Derived3>
    Return_t df_sym(const Eigen::MatrixBase<Derived0> & cos_thetas,
                    const Eigen::MatrixBase<Derived1> & distances,
                    const Eigen::MatrixBase<Derived2> & /*cutoff_vals*/,
                    const Eigen::MatrixBase<Derived3> & /*cutoff_derivatives*/,
                    const std::array<size_t, 3> & ordering = {0, 1, 2}) const {
      auto && cos_theta{cos_thetas(ordering[0])};
      auto && lam_cos_theta_1{1. + this->lambda * cos_theta};
      auto && angular_contrib{math::pow(lam_cos_theta_1, this->zeta)};

      auto && r_ij{distances(ordering[0])};
      auto && r_jk{distances(ordering[1])};
      auto && r_ki{distances(ordering[2])};

      auto && r_ij2{r_ij * r_ij};
      auto && r_ki2{r_ki * r_ki};

      auto && exp_contrib{exp(-this->eta * (r_ij2 + r_ki2))};
      auto && fun_val{this->prefactor * angular_contrib * exp_contrib};

      auto && d_dr_ij{-2 * this->eta * r_ij * fun_val +
                      this->zeta * (this->lambda / r_ki - cos_theta / r_ij) *
                          fun_val / lam_cos_theta_1};
      auto && d_dr_ki{-2 * this->eta * r_ki * fun_val +
                      this->zeta * (this->lambda / r_ij - cos_theta / r_ki) *
                          fun_val / lam_cos_theta_1};
      auto && d_dr_jk{-this->lambda * r_jk * this->zeta * fun_val /
                      (r_ij * r_ki * lam_cos_theta_1)};
      return Return_t(fun_val, {d_dr_ij, d_dr_jk, d_dr_ki});
    }

   protected:
    const json params;
    double zeta;
    double lambda;
    double eta;
    double prefactor;  // prefactor (2^(1-ζ))
  };

  constexpr size_t SymmetryFunction<SymmetryFunctionType::AngularWide>::Order;

  /**
   * there are three angles in triplets, hence the magic number 3
   */
  template <class Manager>
  const Property<double, TripletOrder, Manager, 3> &
  get_cos_angles(Manager & manager) {
    // 1) get or create
    const std::string identifier{"cosines_of_angles"};
    constexpr bool Validate{false}, AllowCreation{true};

    // there are three angles in triplets
    auto & property{
        *manager
             .template get_property<Property<double, TripletOrder, Manager, 3>>(
                 identifier, Validate, AllowCreation)};

    // 2) fresh?
    //   a) all done
    if (property.is_updated()) {
      return property;
    }
    //   b) compute
    auto & triplet_direction_vectors{manager.get_triplet_direction_vectors()};
    for (auto && atom : manager) {
      for (auto && trip : atom.triplets()) {
        auto && vectors{triplet_direction_vectors[trip]};
        auto && v_ij{vectors.col(0)};
        auto && v_jk{vectors.col(1)};
        auto && v_ki{vectors.col(2)};

        property[trip](0) = v_ij.dot(-v_ki);
        property[trip](1) = v_jk.dot(-v_ij);
        property[trip](2) = v_ki.dot(-v_jk);
      }
    }

    property.set_updated_status(true);
    return property;
  }
}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_SYMMETRY_FUNCTIONS_HH_
