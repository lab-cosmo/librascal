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

  using interpolator_fixtures =
      boost::mpl::list<InterpolatorFixture<IntpScalarUniformCubicSpline>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error;
    std::function<double(double)> func;
    std::string function_name;

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
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

      // Checks if interpolator satisfies the given error bound. This condition
      // should be always fulfilled.
      error = compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, Fix::error_bound);
    }
  }

  // Because the hyp1f1 and radial contribution has to be initiated, these tests
  // cannot be included in the function map and have to be specified in a
  // separate tests
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(hyp1f1_interpolator_tests, Fix,
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

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_derivative_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error, intp_error, finite_diff_error;
    std::function<double(double)> func;
    std::function<double(double)> derivative_func;
    std::string function_name;
    std::shared_ptr<typename Fix::Interpolator_t> intp{};

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
      derivative_func = Fix::derivatives[function_name];
      if (verbose) {
        std::cout << "Test for function: " << function_name << std::endl;
      }
      for (bool clamped_boundary_condition : {true, false}) {
        if (verbose) {
          std::cout << "Interpolator interpolator clamped_boundary_condition : "
                    << clamped_boundary_condition << std::endl;
        }
        if (clamped_boundary_condition) {
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
        error = compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
        BOOST_CHECK_LE(error, Fix::error_bound);

        // Checks if using the interpolators derivative intp' is at least as
        // accurate as using finite methods with the interpolators function intp
        // or if they are close. This condition should be always fulfilled.

        // error of using intp'(x_i)
        intp_error = compute_intp_derivative_error<AbsoluteErrorMethod>(
            intp, derivative_func, ref_points);
        // error of using ( intp(x_{i+1})-intp(x_i) ) / ( x_{i+1} - x_i )
        finite_diff_error = compute_intp_finite_diff_error<AbsoluteErrorMethod>(
            intp, derivative_func, ref_points);
        BOOST_CHECK((finite_diff_error - intp_error) > -tol);
        if (verbose) {
          std::cout << "Interpolator derivative compared to interpolator "
                       "function finite method error on test grid: "
                    << finite_diff_error - intp_error << std::endl;
        }

        // The following conditions do not necessary have to be fulfilled, but
        // an error here could point out a bug in the interpolator or question
        // the sanity of it. For an unfulfilled condition the error should be
        // at least be comparatively small.

        // Checks if derivative holds for a slightly higher error bound. Since
        // the order of the error does not increase, it should be close to the
        // error bound. A rigorous error bound can be found in
        // https://doi.org/10.1016/0021-9045(82)90041-7 and could be
        // implemented, but this is good enough.
        error = compute_intp_derivative_error<AbsoluteErrorMethod>(
            intp, derivative_func, ref_points);
        BOOST_CHECK_LE(error, 6 * Fix::error_bound);

        // Check if derivative of the interpolator is as precise as the finite
        // difference method for the points used by the interpolator.
        intp_error = compute_intp_derivative_error<AbsoluteErrorMethod>(
            intp, derivative_func, intp->get_grid_ref());
        finite_diff_error = compute_finite_diff_error<AbsoluteErrorMethod>(
            func, derivative_func, intp->get_grid_ref());
        BOOST_CHECK((finite_diff_error - intp_error) >= 0 ||
                    std::abs(finite_diff_error - intp_error) < tol);
        if (verbose) {
          std::cout << "Interpolator derivative compared to function finite "
                       "method error on interpolator grid: "
                    << std::abs(finite_diff_error - intp_error) << std::endl;
        }
      }
    }
  }

  using IntpScalarUniformCubicSplineRelativeError =
      math::InterpolatorScalarUniformCubicSpline<
          math::RefinementMethod_t::Exponential,
          math::ErrorMethod<math::ErrorMetric_t::Relative>>;
  using interpolator_relative_error_fixtures = boost::mpl::list<
      InterpolatorFixture<IntpScalarUniformCubicSplineRelativeError>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_interpolator_relative_error_tests,
                                   Fix, interpolator_relative_error_fixtures,
                                   Fix) {
    bool verbose{false};

    Vector_t ref_points =
        Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error;
    std::function<double(double)> func;
    std::string function_name;

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
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

      // Checks if interpolator satisfies the given error bound. This condition
      // should be always achieved. The error_bound for the exp function is not
      // completely achieved therefore the 2* factor. Since the intepolator only
      // estimates an error this is still within a reasonable result.
      error = compute_intp_error<RelativeErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, 2 * Fix::error_bound);
    }
  }

  using IntpMatrixUniformCubicSpline =
      math::InterpolatorMatrixUniformCubicSpline<
          math::RefinementMethod_t::Exponential>;

  // This test compares the vectorized interpolator with the scalar
  // interpolators results for the functions map. The results should match
  BOOST_FIXTURE_TEST_CASE(functions_interpolator_vectorized_test,
                          InterpolatorFixture<IntpMatrixUniformCubicSpline>) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2);
    double error;
    std::function<double(double)> func;
    std::function<Matrix_t(double)> vector_func;
    std::string function_name;

    for (const auto pair : functions) {
      function_name = pair.first;
      if (verbose) {
        std::cout << "Test for function: " << function_name << std::endl;
      }
      func = functions[function_name];
      vector_func = [&](double x) {
        Matrix_t mat(1, 1);
        mat << func(x);
        return mat;
      };

      Matrix_t tmp_mat = vector_func(x1);
      int cols{static_cast<int>(tmp_mat.cols())};
      int rows{static_cast<int>(tmp_mat.rows())};
      auto intp{std::make_shared<IntpMatrixUniformCubicSpline>(
          vector_func, x1, x2, error_bound, cols, rows)};
      if (verbose) {
        std::cout << "Interpolator interpolator fineness: "
                  << intp->get_degree_of_fineness() << std::endl;
      }
      if (verbose) {
        std::cout << "Interpolator grid size: " << intp->get_grid_size()
                  << std::endl;
      }

      // Checks if interpolator satisfies the given error bound. This condition
      // should be always fulfilled.
      error = compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, error_bound);
    }
  }

  BOOST_FIXTURE_TEST_CASE(functions_dervative_interpolator_vectorized_test,
                          InterpolatorFixture<IntpMatrixUniformCubicSpline>) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2);
    double error, intp_error, finite_diff_error;
    std::function<double(double)> func;
    std::function<Matrix_t(double)> vector_func;
    std::function<double(double)> derivative_func;
    std::string function_name;

    for (const auto pair : functions) {
      function_name = pair.first;
      if (verbose) {
        std::cout << "Test for function: " << function_name << std::endl;
      }
      func = functions[function_name];
      vector_func = [&](double x) {
        Matrix_t mat(1, 1);
        mat << func(x);
        return mat;
      };
      derivative_func = derivatives[function_name];

      Matrix_t tmp_mat = vector_func(x1);
      int cols = tmp_mat.cols();
      int rows = tmp_mat.rows();
      auto intp{std::make_shared<IntpMatrixUniformCubicSpline>(
          vector_func, x1, x2, error_bound, cols, rows)};
      if (verbose) {
        std::cout << "Interpolator interpolator fineness: "
                  << intp->get_degree_of_fineness() << std::endl;
      }
      if (verbose) {
        std::cout << "Interpolator grid size: " << intp->get_grid_size()
                  << std::endl;
      }

      // Checks if interpolator satisfies the given error bound. This condition
      // should be always fulfilled.
      error = compute_intp_error<AbsoluteErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, error_bound);

      // The following conditions do not necessary have to be fulfilled, but an
      // error here could point out a bug in the interpolator or question the
      // sanity of it or be because the derivative function does not agree with
      // the function. For an unfulfillde condition the error should be at least
      // be comparatively small.

      // Checks if using the interpolators derivative is at least as accurate as
      // using finite methods with the interpolators function or if they are
      // close. This condition should be always fulfilled.
      intp_error = compute_intp_derivative_error<AbsoluteErrorMethod>(
          intp, derivative_func, ref_points);
      finite_diff_error = compute_intp_finite_diff_error<AbsoluteErrorMethod>(
          intp, derivative_func, ref_points);
      BOOST_CHECK((finite_diff_error - intp_error) >= 0 ||
                  std::abs(finite_diff_error - intp_error) < tol);
      if (verbose) {
        std::cout << "Interpolator derivative compared to interpolator "
                     "function finite method error on test grid: "
                  << std::abs(finite_diff_error - intp_error) << std::endl;
      }

      // Checks if derivative holds for a slightly higher error bound. Since the
      // order of the error does not increase, it should be close to the error
      // bound. A rigorous error bound can be found in
      // https://doi.org/10.1016/0021-9045(82)90041-7 and could be implemented,
      // but this is good enough.
      error = compute_intp_derivative_error<AbsoluteErrorMethod>(
          intp, derivative_func, ref_points);
      BOOST_CHECK_LE(error, 6 * error_bound);

      // Check if derivative of the interpolator is as precise as the finite
      // difference method for the points used by the interpolator.
      intp_error = compute_intp_derivative_error<AbsoluteErrorMethod>(
          intp, derivative_func, intp->get_grid_ref());
      finite_diff_error = compute_finite_diff_error<AbsoluteErrorMethod>(
          func, derivative_func, intp->get_grid_ref());
      BOOST_CHECK((finite_diff_error - intp_error) >= 0 ||
                  std::abs(finite_diff_error - intp_error) < tol);
      if (verbose) {
        std::cout << "Interpolator derivative compared to function finite "
                     "method error on interpolator grid: "
                  << std::abs(finite_diff_error - intp_error) << std::endl;
      }
    }
  }

  // Because the radial contribution has to be initiated, this test is not
  // included in the function map and has to be specified in a separate test. It
  // checks if the standard interpolator and vectorized interpolator return the
  // same results and if the vectorized interpolator is fulfills its error bound
  BOOST_FIXTURE_TEST_CASE(radial_contribution_interpolator_vectorized_test,
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
