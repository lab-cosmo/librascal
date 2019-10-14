/**
 * file   test_math_interpolator.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   17 June 2019
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

#include "test_math_interpolator.hh"

namespace rascal {

  // TODO(all) not sure about the naming convention of tests, camelcase ...
  BOOST_AUTO_TEST_SUITE(MathInterpolatorTests);

  using IntpScalarUniformCubicSpline =
      math::InterpolatorScalarUniformCubicSpline<
          math::RefinementMethod_t::Exponential>;

  using IntpScalarUniformCubicSplineRelativeError =
      math::InterpolatorScalarUniformCubicSpline<
          math::RefinementMethod_t::Exponential,
          math::ErrorMethod<math::ErrorMetric_t::Relative>>;

  using IntpMatrixUniformCubicSpline =
      math::InterpolatorMatrixUniformCubicSpline<
          math::RefinementMethod_t::Exponential>;

  using interpolator_fixtures = boost::mpl::list<
      InterpolatorFixture<IntpScalarUniformCubicSpline>,
      InterpolatorFixture<IntpScalarUniformCubicSplineRelativeError>,
      InterpolatorFixture<IntpMatrixUniformCubicSpline>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(interpolator_constructor_test, Fix,
                                   interpolator_fixtures, Fix) {
    auto intp{std::make_shared<typename Fix::Interpolator_t>(
        Fix::functions["identity"], Fix::x1, Fix::x2, Fix::error_bound)};
  }

  /**
   * Tests for scalar functions.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_interpolator_test, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error;
    std::function<double(double)> func;
    std::string function_name;

    for (const auto & it : Fix::functions) {
      function_name = it.first;
      func = it.second;
      if (verbose) {
        std::cout << "Test for function: " << function_name << std::endl;
      }

      auto intp{std::make_shared<typename Fix::Interpolator_t>(
          func, Fix::x1, Fix::x2, Fix::error_bound)};
      if (verbose) {
        std::cout << "Interpolator interpolator fineness: "
                  << intp->get_degree_of_fineness() << std::endl;
      }
      if (verbose) {
        std::cout << "Interpolator grid size: " << intp->get_grid_size()
                  << std::endl;
      }

      // Checks if interpolator satisfies the given error bound. The error_bound
      // for the exp function is not completely achieved using a RelativeError
      // function therefore the 2* factor. Since the intepolator only estimates
      // an error this is still within a reasonable result.
      error = compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, 2 * Fix::error_bound);
    }
  }

  /**
   * Because the hyp1f1 and radial contribution has to be initialized, this
   * test cannot be included in `functions_interpolator_test`, therefore these
   * two get separate tests.
   *
   * Hyp1f1
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(hyp1f1_interpolator_test, Fix,
                                   interpolator_fixtures, Fix) {
    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    std::function<double(double)> func = [&](double x) {
      return Fix::hyp1f1.calc(x);
    };
    auto intp{std::make_shared<typename Fix::Interpolator_t>(
        func, Fix::x1, Fix::x2, Fix::error_bound)};
    double error =
        compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
    BOOST_CHECK_LE(error, Fix::error_bound);
  }

  /*
   * RadialContribution
   */
  BOOST_FIXTURE_TEST_CASE(radial_contribution_test,
                          InterpolatorFixture<IntpScalarUniformCubicSpline>) {
    Vector_t ref_points =
        Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    std::function<double(double)> func = [&](double x) {
      return this->radial_contr
          .compute_contribution<rascal::internal::AtomicSmearingType::Constant>(
              x, 0.5)(0, 0);
    };
    auto intp{std::make_shared<IntpScalarUniformCubicSpline>(
        func, this->x1, this->x2, this->error_bound)};
    double error =
        compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
    BOOST_CHECK_LE(error, this->error_bound);
  }

  /**
   * Tests for scalar derivative functions.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_derivative_interpolator_test, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double intp_error, intp_derivative_error, func_finite_diff_error,
        intp_finite_diff_error;
    std::function<double(double)> func;
    std::function<double(double)> derivative_func;
    std::string function_name;
    std::shared_ptr<typename Fix::Interpolator_t> intp{};

    for (const auto & it : Fix::functions) {
      function_name = it.first;
      func = it.second;
      derivative_func = Fix::derivatives[function_name];
      if (verbose) {
        std::cout << "Test for function: " << function_name << std::endl;
      }
      if (verbose) {
        std::cout << "Interpolator interpolator clamped_boundary_conditions : "
                  << Fix::clamped_boundary_conditions << std::endl;
      }
      if (Fix::clamped_boundary_conditions) {
        intp = std::make_shared<typename Fix::Interpolator_t>(
            func, Fix::x1, Fix::x2, Fix::error_bound, 100000, 5, true,
            derivative_func(Fix::x1), derivative_func(Fix::x2));
      } else {
        intp = std::make_shared<typename Fix::Interpolator_t>(
            func, Fix::x1, Fix::x2, Fix::error_bound);
      }
      if (verbose) {
        std::cout << "Interpolator interpolator fineness: "
                  << intp->get_degree_of_fineness() << std::endl;
      }
      if (verbose) {
        std::cout << "Interpolator grid size: " << intp->get_grid_size()
                  << std::endl;
      }

      // Checks if interpolator satisfies the given error bound. This
      // condition should be always fulfilled.
      intp_error =
          compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(intp_error, Fix::error_bound);

      // Checks if using the interpolators derivative intp' is at least as
      // accurate as using finite methods with the interpolators function intp
      // or if they are close. This condition should be always fulfilled.

      // error using intp'(x_i)
      intp_derivative_error =
          compute_intp_derivative_error<AbsoluteErrorMethod>(
              intp, derivative_func, ref_points);
      // error using ( intp(x_{i+1})-intp(x_i) ) / ( x_{i+1} - x_i )
      intp_finite_diff_error =
          compute_intp_finite_diff_error<AbsoluteErrorMethod>(
              intp, derivative_func, ref_points);
      BOOST_CHECK((intp_finite_diff_error - intp_derivative_error) > -tol);
      if (verbose) {
        std::cout << "Interpolator derivative compared to estimation of "
                     "derivative using finite difference with the interpolated"
                     "function on test grid: "
                  << intp_finite_diff_error - intp_derivative_error
                  << std::endl;
      }

      // Checks if derivative holds for a slightly higher error bound. Since
      // the order of the error does not increase, it should be close to the
      // error bound. A rigorous error bound for CubicSpline can be found in
      // https://doi.org/10.1016/0021-9045(82)90041-7 and could be
      // implemented, but this is good enough.
      BOOST_CHECK_LE(intp_derivative_error, 6 * Fix::error_bound);

      // Check if derivative of the interpolator is as least as precise as the
      // finite difference method with the correct function on the grid used
      // by the interpolator.

      // error using intp'(x_i)
      intp_derivative_error =
          compute_intp_derivative_error<AbsoluteErrorMethod>(
              intp, derivative_func, intp->get_grid_ref());
      // error using ( f(x_{i+1})-f(x_i) ) / ( x_{i+1} - x_i )
      func_finite_diff_error = compute_finite_diff_error<AbsoluteErrorMethod>(
          func, derivative_func, intp->get_grid_ref());
      BOOST_CHECK((func_finite_diff_error - intp_derivative_error) > -tol);
      if (verbose) {
        std::cout << "Interpolator derivative compared to estimation of "
                     "the derivative using finite difference with the "
                     "correct function on interpolator grid: "
                  << func_finite_diff_error - intp_derivative_error
                  << std::endl;
      }
    }
  }

  /**
   * The radial contribution test for the matrix interpolator
   */
  BOOST_FIXTURE_TEST_CASE(radial_contribution_matrix_interpolator_test,
                          InterpolatorFixture<IntpMatrixUniformCubicSpline>) {
    Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2);

    std::function<Matrix_t(double)> func = [&](double x) {
      return this->radial_contr.compute_neighbour_contribution(x, 0.5);
    };
    Matrix_t tmp_mat = func(x1);
    int cols = tmp_mat.cols();
    int rows = tmp_mat.rows();
    auto intp{std::make_shared<IntpMatrixUniformCubicSpline>(
        func, x1, x2, error_bound, cols, rows)};

    int matrix_size = max_radial * (max_angular + 1);
    Matrix_t intp_val = Matrix_t::Zero(ref_points.size(), matrix_size);
    Matrix_t intp_ref = Matrix_t::Zero(ref_points.size(), matrix_size);
    for (int i{0}; i < ref_points.size() - 1; i++) {
      intp_val.row(i) = Eigen::Map<Vector_t>(
          intp->interpolate(ref_points(i)).data(), matrix_size);
      intp_ref.row(i) =
          Eigen::Map<Vector_t>(func(ref_points(i)).data(), matrix_size);
    }
    double error{
        (intp_val - intp_ref).array().abs().colwise().mean().maxCoeff()};
    BOOST_CHECK_LE(error, error_bound);
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
