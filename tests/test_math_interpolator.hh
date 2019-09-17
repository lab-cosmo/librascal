/**
 * file   test_math_interpolator.hh
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

#include "tests.hh"
#include "test_math.hh"
#include "math/interpolator.hh"

#ifndef TESTS_TEST_MATH_INTERPOLATOR_HH_
#define TESTS_TEST_MATH_INTERPOLATOR_HH_

namespace rascal {
  using Vector_t = math::Vector_t;
  using Matrix_t = math::Matrix_t;
  using Vector_Ref = math::Vector_Ref;
  using Matrix_Ref = math::Matrix_Ref;

  using AbsoluteMeanErrorMethod =
      math::ErrorMethod<math::ErrorMetric_t::Absolute, math::ErrorNorm_t::Mean>;
  using RelativeMeanErrorMethod =
      math::ErrorMethod<math::ErrorMetric_t::Relative, math::ErrorNorm_t::Mean>;

  static double convert_to_scalar(const double x) { return x; }
  static double convert_to_scalar(const math::Matrix_t & mat) {
    return mat(0, 0);
  }

  template <class Interpolator>
  static double intp_ref_mean_error(std::shared_ptr<Interpolator> & intp,
                                    const math::Vector_Ref & ref_points) {
    return (intp->interpolate(ref_points) - intp->eval(ref_points))
        .array()
        .abs()
        .mean();
  }
  // These functions calculate the error between the interpolated function and
  // the interpolation on a set of reference points for different

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

  template <class ErrorMethod, class Interpolator>
  static double
  compute_intp_derivative_error(std::shared_ptr<Interpolator> & intp,
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

  template <class ErrorMethod, class Interpolator>
  static double
  compute_intp_finite_diff_error(std::shared_ptr<Interpolator> & intp,
                                 std::function<double(double)> derivative_func,
                                 const math::Vector_Ref & ref_points) {
    Vector_t intp_refs = Vector_t::Zero(ref_points.size() - 1);
    Vector_t finite_diff_derivative = Vector_t::Zero(ref_points.size() - 1);
    for (int i{0}; i < ref_points.size() - 2; i++) {
      intp_refs(i) = derivative_func(ref_points(i));
      finite_diff_derivative(i) =
          (convert_to_scalar(intp->interpolate(ref_points(i + 1))) -
           convert_to_scalar(intp->interpolate(ref_points(i)))) /
          (ref_points(i + 1) - ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(finite_diff_derivative),
                                             Vector_Ref(intp_refs));
  }

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

  template <class Interpolator>
  struct InterpolatorFixture {
    InterpolatorFixture<Interpolator>() {}
    typedef Interpolator Interpolator_t;

    double x1{0};
    double x2{1};
    double error_bound{1e-5};
    int nb_ref_points{100};

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

    math::Hyp1f1 hyp1f1{0.5 * (20 + 20 + 3), 20 + 1.5, 200, 1e-15};
    const int max_radial{3};
    const int max_angular{max_radial - 1};
    const json fc_hypers{{"type", "Constant"},
                         {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}};
    const json hypers{
        {"gaussian_density", fc_hypers},
        {"max_radial", max_radial},
        {"max_angular", max_angular},
        {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};

    // a=0.5*(n+l+3), b = l+1.5, mmax, tolerance
    internal::RadialContribution<rascal::internal::RadialBasisType::GTO>
        radial_contr{internal::RadialContribution<
            rascal::internal::RadialBasisType::GTO>(hypers)};
  };

}  // namespace rascal

#endif  // TESTS_TEST_MATH_INTERPOLATOR_HH_
