/**
 * @file   test_math_interpolator.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   26 August 2019
 *
 * @brief Tests the implementation of interpolator
 *
 * Copyright © 2019  Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "math/interpolator.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "test_math.hh"

#include <boost/mpl/list.hpp>

#ifndef TESTS_TEST_MATH_INTERPOLATOR_HH_
#define TESTS_TEST_MATH_INTERPOLATOR_HH_

namespace rascal {
  using Vector_t = math::Vector_t;
  using Matrix_t = math::Matrix_t;
  using Vector_Ref = math::Vector_Ref;
  using Matrix_Ref = math::Matrix_Ref;

  using AbsoluteErrorMethod = math::ErrorMethod<math::ErrorMetric_t::Absolute>;
  using RelativeErrorMethod = math::ErrorMethod<math::ErrorMetric_t::Relative>;

  /**
   * For the `functions_matrix_interpolator_test` test we can reuse the
   * error functions for the scalar interpolator by converting the matrix of
   * shape (1,1) to a scalar value before.
   */
  static double convert_to_scalar(const double x) { return x; }
  static double convert_to_scalar(const math::Matrix_t & mat) {
    return mat(0, 0);
  }

  /**
   * This function calculate the error between the interpolated function
   * of `intp` and the correct function `ref_func` on the points `ref_points`.
   *
   * @paramt ErrorMethod the error method calculating reducing grid error to a
   *         scalar value.
   * @param intp scalar interpolator or matrix interpolator for a matrix of
   *        shape (1,1).
   * @param ref_func correct function.
   * @param ref_points reference grid points to estimate the error.
   */
  template <class ErrorMethod, class Interpolator>
  static double compute_intp_error(std::shared_ptr<Interpolator> & intp,
                                   std::function<double(double)> ref_func,
                                   const math::Vector_Ref & ref_points) {
    math::Vector_t intp_vals = math::Vector_t::Zero(ref_points.size());
    math::Vector_t intp_refs = math::Vector_t::Zero(ref_points.size());
    for (int i{0}; i < ref_points.size() - 1; i++) {
      intp_vals(i) = convert_to_scalar(intp->interpolate(ref_points(i)));
      intp_refs(i) = ref_func(ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(intp_vals),
                                             Vector_Ref(intp_refs));
  }

  /**
   * This function calculate the error between the derivative of the
   * interpolated function of `intp` and the correct derivative
   * `derivative_func` on the points `ref_points`.
   *
   * @paramt ErrorMethod the error method calculating reducing grid error to a
   *         scalar value.
   * @param intp scalar interpolator or matrix interpolator for a matrix of
   *        shape (1,1).
   * @param derivative_func correct derivative.
   * @param ref_points reference grid points to estimate the error.
   * @return error
   */
  template <class ErrorMethod, class Interpolator>
  static double
  compute_intp_derivative_error(const std::shared_ptr<Interpolator> & intp,
                                std::function<double(double)> derivative_func,
                                const math::Vector_Ref & ref_points) {
    math::Vector_t intp_vals = math::Vector_t::Zero(ref_points.size());
    math::Vector_t intp_refs = math::Vector_t::Zero(ref_points.size());
    for (int i{0}; i < ref_points.size() - 1; i++) {
      intp_vals(i) =
          convert_to_scalar(intp->interpolate_derivative(ref_points(i)));
      intp_refs(i) = derivative_func(ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(intp_vals),
                                             Vector_Ref(intp_refs));
  }

  /**
   * This function calculate the error between the derivative estimation of the
   * interpolated function of `intp` using finite difference method and the
   * correct derivative `derivative_func` on the points `ref_points`.
   *
   * @paramt ErrorMethod the error method calculating reducing grid error to a
   *         scalar value.
   * @param intp scalar interpolator or matrix interpolator for a matrix of
   *        shape (1,1).
   * @param derivative_func correct derivative.
   * @param ref_points reference grid points to estimate the error.
   */
  template <class ErrorMethod, class Interpolator>
  static double
  compute_intp_finite_diff_error(const std::shared_ptr<Interpolator> & intp,
                                 std::function<double(double)> derivative_func,
                                 const math::Vector_Ref & ref_points) {
    Vector_t derivative_refs = Vector_t::Zero(ref_points.size() - 1);
    Vector_t intp_finite_diff_derivative =
        Vector_t::Zero(ref_points.size() - 1);
    for (int i{0}; i < ref_points.size() - 2; i++) {
      derivative_refs(i) = derivative_func(ref_points(i));
      // finite difference method
      intp_finite_diff_derivative(i) =
          (convert_to_scalar(intp->interpolate(ref_points(i + 1))) -
           convert_to_scalar(intp->interpolate(ref_points(i)))) /
          (ref_points(i + 1) - ref_points(i));
    }
    return ErrorMethod::compute_global_error(
        Vector_Ref(intp_finite_diff_derivative), Vector_Ref(derivative_refs));
  }

  /**
   * This function calculate the error between the derivative estimation of the
   * function `func` using finite difference method and the
   * correct derivative `derivative_func`.
   *
   * @paramt ErrorMethod the error method calculating the grid error and
   *         reducing it to a scalar value.
   * @param func the correct function
   * @param derivative_func correct derivative.
   * @param ref_points reference grid points to estimate the error.
   * @return error
   */
  template <class ErrorMethod>
  static double
  compute_finite_diff_error(std::function<double(double)> func,
                            std::function<double(double)> derivative_func,
                            const math::Vector_Ref & ref_points) {
    Vector_t derivative_refs = Vector_t::Zero(ref_points.size() - 1);
    Vector_t finite_diff_derivative = Vector_t::Zero(ref_points.size() - 1);
    for (int i{0}; i < ref_points.size() - 2; i++) {
      derivative_refs(i) = derivative_func(ref_points(i));
      finite_diff_derivative(i) =
          (func(ref_points(i + 1)) - func(ref_points(i))) /
          (ref_points(i + 1) - ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(finite_diff_derivative),
                                             Vector_Ref(derivative_refs));
  }

  template <class Interpolator, bool ClampedBoundaryConditions = false>
  struct InterpolatorFixture {
    using RadialIntegral_t =
        internal::RadialContribution<internal::RadialBasisType::GTO>;
    InterpolatorFixture<Interpolator, ClampedBoundaryConditions>() {}
    typedef Interpolator Interpolator_t;
    static constexpr bool clamped_boundary_conditions{
        ClampedBoundaryConditions};

    double x1{0};
    double x2{1};
    double error_bound{1e-5};
    int nb_ref_points{100};

    // different functions within the interval [-1,1] for the range [x1,x2] are
    // used for testing, because the absolute error is mainly used
    std::map<std::string, std::function<double(double)>> functions{
        {"identity", {[](double x) { return x; }}},
        {"polynomial", {[](double x) { return x * x * x / 900 - x / 2 + 3; }}},
        {"exp", {[](double x) { return std::exp(0.5 * x); }}},
        {"sin", {[](double x) { return std::sin(x); }}}};
    std::map<std::string, std::function<double(double)>> derivatives{
        {"identity", {[](double) { return 1; }}},
        {"polynomial", {[](double x) { return x * x / 300 - 0.5; }}},
        {"exp", {[](double x) { return 0.5 * std::exp(0.5 * x); }}},
        {"sin", {[](double x) { return std::cos(x); }}}};

    // For a == b the hy1f1 function is the same as exp therefore we test the
    // case a != b radial = 5; angular = 4; a = 0.5 * (radial + angular + 3); b
    // = angular + 1.5; math::Hyp1f1 hyp1f1{a, b, 200, 1e-15};
    math::Hyp1f1 hyp1f1{0.5 * (5 + 4 + 3), 4 + 1.5, 200, 1e-15};
    const int max_radial{3};
    const int max_angular{max_radial - 1};
    const json fc_hypers{{"type", "Constant"},
                         {"gaussian_sigma", {{"value", 0.5}, {"unit", "AA"}}}};
    const json hypers{
        {"gaussian_density", fc_hypers},
        {"max_radial", max_radial},
        {"max_angular", max_angular},
        {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "AA"}}}}}};

    // a=0.5*(n+l+3), b = l+1.5, mmax, tolerance
    RadialIntegral_t radial_contr{RadialIntegral_t(hypers)};
  };

}  // namespace rascal

#endif  // TESTS_TEST_MATH_INTERPOLATOR_HH_
