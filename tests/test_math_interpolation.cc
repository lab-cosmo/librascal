/**
 * file   test_math_interpolation.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   17 June 2019
 *
 * @brief Test the implementation of Interpolator 
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

namespace rascal {
  namespace math {
  using AbsoluteMeanErrorMethod = ErrorMethod<ErrorMetric_t::Absolute, ErrorNorm_t::Mean>;
  using RelativeMeanErrorMethod = ErrorMethod<ErrorMetric_t::Relative, ErrorNorm_t::Mean>;

  static double convert_to_scalar(const double x) {
    return x;
  }
  static double convert_to_scalar(const Matrix_t & mat) {
    return mat(0,0);
  }

  template <class Interpolator>
  static double intp_ref_mean_error(Interpolator & intp, const math::Vector_Ref & ref_points) {
    return (intp.interpolate(ref_points) - intp.eval(ref_points)).array().abs().mean();
  }
  // These functions calculate the error between the interpolated function and the interpolation on a set of reference points for different 

  template <class ErrorMethod, class Interpolator>
  static double compute_intp_error(Interpolator & intp, std::function<double(double)> ref_func, const math::Vector_Ref & ref_points) {
    math::Vector_t intp_vals = math::Vector_t::Zero(ref_points.size());
    math::Vector_t intp_refs = math::Vector_t::Zero(ref_points.size());
    for (int i{0}; i<ref_points.size()-1; i++) {
      intp_vals(i) = convert_to_scalar(intp.interpolate(ref_points(i)));
      intp_refs(i) = ref_func(ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(intp_vals), Vector_Ref(intp_refs));
  }


  template <class ErrorMethod, class Interpolator>
  static double compute_intp_derivative_error(Interpolator & intp, std::function<double(double)> derivative_func, const math::Vector_Ref & ref_points) {
    math::Vector_t intp_vals = math::Vector_t::Zero(ref_points.size());
    math::Vector_t intp_refs = math::Vector_t::Zero(ref_points.size());
    for (int i{0}; i<ref_points.size()-1; i++) {
      intp_vals(i) = convert_to_scalar(intp.interpolate_derivative(ref_points(i)));
      intp_refs(i) = derivative_func(ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(intp_vals), Vector_Ref(intp_refs));
  }


  template <class ErrorMethod, class Interpolator>
  static double compute_intp_finite_diff_error(Interpolator & intp, std::function<double(double)> derivative_func, const math::Vector_Ref & ref_points) {
    Vector_t intp_refs = Vector_t::Zero(ref_points.size()-1);
    Vector_t finite_diff_derivative = Vector_t::Zero(ref_points.size()-1);
    for (int i{0}; i<ref_points.size()-2; i++) {
      intp_refs(i) = derivative_func(ref_points(i));
      finite_diff_derivative(i) = (convert_to_scalar(intp.interpolate(ref_points(i+1)))-convert_to_scalar(intp.interpolate(ref_points(i))))/(ref_points(i+1)-ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(finite_diff_derivative), Vector_Ref(intp_refs));
  }

  template <class ErrorMethod>
  static double compute_finite_diff_error(std::function<double(double)> func, std::function<double(double)> derivative_func, const math::Vector_Ref & ref_points) {
    Vector_t derivative_refs = Vector_t::Zero(ref_points.size()-1);
    Vector_t finite_diff_derivative = Vector_t::Zero(ref_points.size()-1);
    for (int i{0}; i<ref_points.size()-2; i++) {
      derivative_refs(i) = derivative_func(ref_points(i));
      finite_diff_derivative(i) = (func(ref_points(i+1))-func(ref_points(i)))/(ref_points(i+1)-ref_points(i));
    }
    return ErrorMethod::compute_global_error(Vector_Ref(finite_diff_derivative), Vector_Ref(derivative_refs));
  }
  
  //
  // ----------------------------------------------
  //

  BOOST_AUTO_TEST_SUITE(MathInterpolatorTests);

  using InterpolatorDefault = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;

  using interpolator_fixtures = boost::mpl::list<
                   InterpolatorFixture<InterpolatorDefault>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error;
    std::function<double(double)> func;
    std::string function_name;

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
      if (verbose) { std::cout << "Test for function: " << function_name << std::endl;}

      Fix::intp.initialize(func, Fix::x1, Fix::x2, Fix::error_bound);
      if (verbose) { std::cout << "Interpolator interpolator fineness: " << Fix::intp.fineness << std::endl;}
      if (verbose) { std::cout << "Interpolator grid size: " << Fix::intp.grid.size() << std::endl;}

      // Checks if interpolator satisfies the given error bound. This condition should be always fulfilled.
      error = compute_intp_error<AbsoluteMeanErrorMethod>(Fix::intp, func, ref_points);
      BOOST_CHECK_LE(error, Fix::error_bound);
    }
  }

  // Because the hyp1f1 and radial contribution has to be initiated, these tests cannot be included in the function map and have to be specified in a separate tests
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(hyp1f1_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    Vector_t ref_points = Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    std::function<double(double)> func  = [&](double x) {return Fix::hyp1f1.calc(x);};
    Fix::intp.initialize(func, Fix::x1, Fix::x2, Fix::error_bound);
    double error = compute_intp_error<AbsoluteMeanErrorMethod>(Fix::intp, func, ref_points);
    BOOST_CHECK_LE(error, Fix::error_bound);
  }
  
  BOOST_FIXTURE_TEST_CASE(radial_contribution_test, InterpolatorFixture<InterpolatorDefault>) {
    Vector_t ref_points = Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    std::function<double(double)> func =
        [&](double x) {return this->radial_contr.compute_contribution<rascal::internal::AtomicSmearingType::Constant>(x, 0.5)(0,0);};
    this->intp.initialize(func, this->x1, this->x2, this->error_bound);
    double error = compute_intp_error<AbsoluteMeanErrorMethod>(this->intp, func, ref_points);
    BOOST_CHECK_LE(error, this->error_bound);
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_derivative_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error, intp_error, finite_diff_error;
    std::function<double(double)> func;
    std::function<double(double)> derivative_func;
    std::string function_name;

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
      derivative_func = Fix::derivatives[function_name];
      if (verbose) { std::cout << "Test for function: " << function_name << std::endl;}

      for (bool clamped_boundary_condition : {true, false}) {
        if (verbose) { std::cout << "Interpolator interpolator clamped_boundary_condition : " << clamped_boundary_condition << std::endl;}
        if (clamped_boundary_condition) {
          Fix::intp.initialize(func, Fix::x1, Fix::x2, Fix::error_bound,
              derivative_func(Fix::x1), derivative_func(Fix::x2));
        } else {
          Fix::intp.initialize(func, Fix::x1, Fix::x2, Fix::error_bound);
        }
        if (verbose) { std::cout << "Interpolator interpolator fineness: " << Fix::intp.fineness << std::endl;}
        if (verbose) { std::cout << "Interpolator grid size: " << Fix::intp.grid.size() << std::endl;}

        // Checks if interpolator satisfies the given error bound. This condition should be always fulfilled.
        error = compute_intp_error<AbsoluteMeanErrorMethod>(Fix::intp, func, ref_points);
        BOOST_CHECK_LE(error, Fix::error_bound);

       // The following conditions do not necessary have to be fulfilled, but an error here could point out a bug in the interpolator or question the sanity of it or be because the derivative function does not agree with the function. For an unfulfillde condition the error should be at least be comparatively small.

        // Checks if using the interpolators derivative is at least as accurate as using finite methods with the interpolators function or if they are close. This condition should be always fulfilled.
        intp_error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(Fix::intp, derivative_func, ref_points);
        finite_diff_error = compute_intp_finite_diff_error<AbsoluteMeanErrorMethod>(Fix::intp, derivative_func, ref_points);
        BOOST_CHECK((finite_diff_error - intp_error) >= 0 || std::abs(finite_diff_error - intp_error) < tol);
        if (verbose) { std::cout << "Interpolator derivative compared to interpolator function finite method error on test grid: " << std::abs(finite_diff_error - intp_error) << std::endl;}

        // Checks if derivative holds for a slightly higher error bound. Since the order of the error does not increase, it should be close to the error bound. A rigorous error bound can be found in https://doi.org/10.1016/0021-9045(82)90041-7 and could be implemented, but this is good enough.
        error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(Fix::intp, derivative_func, ref_points);
        BOOST_CHECK_LE(error, 6*Fix::error_bound);

        // Check if derivative of the interpolator is as precise as the finite difference method for the points used by the interpolator. 
        intp_error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(Fix::intp, derivative_func, Fix::intp.grid);
        finite_diff_error = compute_finite_diff_error<AbsoluteMeanErrorMethod>(func, derivative_func, Fix::intp.grid);
        BOOST_CHECK((finite_diff_error - intp_error) >= 0 || std::abs(finite_diff_error - intp_error) < tol);
        if (verbose) { std::cout << "Interpolator derivative compared to function finite method error on interpolator grid: " << std::abs(finite_diff_error - intp_error) << std::endl;}
      }
    }
  }

  using InterpolatorUniformRelativeError = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>,
      ErrorMethod<ErrorMetric_t::Relative, ErrorNorm_t::Mean>
        >;
  using interpolator_relative_error_fixtures = boost::mpl::list<
                   InterpolatorFixture<InterpolatorUniformRelativeError>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(functions_interpolator_relative_error_tests, Fix,
                                   interpolator_relative_error_fixtures , Fix) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2);
    double error;
    std::function<double(double)> func;
    std::string function_name;

    for (const auto pair : Fix::functions) {
      function_name = pair.first;
      func = Fix::functions[function_name];
      if (verbose) { std::cout << "Test for function: " << function_name << std::endl;}

      Fix::intp.initialize(func, Fix::x1, Fix::x2, Fix::error_bound);
      if (verbose) { std::cout << "Interpolator interpolator fineness: " << Fix::intp.fineness << std::endl;}
      if (verbose) { std::cout << "Interpolator grid size: " << Fix::intp.grid.size() << std::endl;}

      // Checks if interpolator satisfies the given error bound. This condition should be always fulfilled.
      error = compute_intp_error<ErrorMethod<ErrorMetric_t::Relative, ErrorNorm_t::Mean>>(Fix::intp, func, ref_points);
      BOOST_CHECK_LE(error, Fix::error_bound);
    }
  }


  using InterpolatorVectorizedDefault = InterpolatorVectorized<
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;

  // This test compares the vectorized interpolator with the scalar interpolators results for the functions map. The results should match
  BOOST_FIXTURE_TEST_CASE(functions_interpolator_vectorized_test,
        InterpolatorFixture<InterpolatorVectorizedDefault>) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2);
    double error;
    std::function<double(double)> func;
    std::function<Matrix_t(double)> vectorized_func;
    std::string function_name;

    for (const auto pair : functions) {
      function_name = pair.first;
      if (verbose) { std::cout << "Test for function: " << function_name << std::endl;}
      func = functions[function_name];
      vectorized_func = [&](double x) {Matrix_t mat(1,1); mat << func(x); return mat;};

      intp.initialize(vectorized_func, x1, x2, error_bound);
      if (verbose) { std::cout << "Interpolator interpolator fineness: " << intp.fineness << std::endl;}
      if (verbose) { std::cout << "Interpolator grid size: " << intp.grid.size() << std::endl;}

      // Checks if interpolator satisfies the given error bound. This condition should be always fulfilled.
      error = compute_intp_error<AbsoluteMeanErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, error_bound);
    }
  }

  BOOST_FIXTURE_TEST_CASE(functions_dervative_interpolator_vectorized_test,
        InterpolatorFixture<InterpolatorVectorizedDefault>) {
    bool verbose{false};

    Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2);
    double error, intp_error, finite_diff_error;
    std::function<double(double)> func;
    std::function<Matrix_t(double)> vectorized_func;
    std::function<double(double)> derivative_func;
    std::string function_name;

    for (const auto pair : functions) {
      function_name = pair.first;
      if (verbose) { std::cout << "Test for function: " << function_name << std::endl;}
      func = functions[function_name];
      vectorized_func = [&](double x) {Matrix_t mat(1,1); mat << func(x); return mat;};
      derivative_func = derivatives[function_name];

      intp.initialize(vectorized_func, x1, x2, error_bound);
      if (verbose) { std::cout << "Interpolator interpolator fineness: " << intp.fineness << std::endl;}
      if (verbose) { std::cout << "Interpolator grid size: " << intp.grid.size() << std::endl;}

      // Checks if interpolator satisfies the given error bound. This condition should be always fulfilled.
      error = compute_intp_error<AbsoluteMeanErrorMethod>(intp, func, ref_points);
      BOOST_CHECK_LE(error, error_bound);

     // The following conditions do not necessary have to be fulfilled, but an error here could point out a bug in the interpolator or question the sanity of it or be because the derivative function does not agree with the function. For an unfulfillde condition the error should be at least be comparatively small.

      // Checks if using the interpolators derivative is at least as accurate as using finite methods with the interpolators function or if they are close. This condition should be always fulfilled.
      intp_error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(intp, derivative_func, ref_points);
      finite_diff_error = compute_intp_finite_diff_error<AbsoluteMeanErrorMethod>(intp, derivative_func, ref_points);
      BOOST_CHECK((finite_diff_error - intp_error) >= 0 || std::abs(finite_diff_error - intp_error) < tol);
      if (verbose) { std::cout << "Interpolator derivative compared to interpolator function finite method error on test grid: " << std::abs(finite_diff_error - intp_error) << std::endl;}

      // Checks if derivative holds for a slightly higher error bound. Since the order of the error does not increase, it should be close to the error bound. A rigorous error bound can be found in https://doi.org/10.1016/0021-9045(82)90041-7 and could be implemented, but this is good enough.
      error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(intp, derivative_func, ref_points);
      BOOST_CHECK_LE(error, 6*error_bound);

      // Check if derivative of the interpolator is as precise as the finite difference method for the points used by the interpolator. 
      intp_error = compute_intp_derivative_error<AbsoluteMeanErrorMethod>(intp, derivative_func, intp.grid);
      finite_diff_error = compute_finite_diff_error<AbsoluteMeanErrorMethod>(func, derivative_func, intp.grid);
      BOOST_CHECK((finite_diff_error - intp_error) >= 0 || std::abs(finite_diff_error - intp_error) < tol);
      if (verbose) { std::cout << "Interpolator derivative compared to function finite method error on interpolator grid: " << std::abs(finite_diff_error - intp_error) << std::endl;}
    }
  }

  // Because the radial contribution has to be initiated, this test is not included in the function map and has to be specified in a separate test. It checks if the standard interpolator and vectorized interpolator return the same results and if the vectorized interpolator is fulfills its error bound
  BOOST_FIXTURE_TEST_CASE(radial_contribution_interpolator_vectorized_test,
        InterpolatorFixture<InterpolatorVectorizedDefault>) {
      Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2); 

      std::function<Matrix_t(double)> func = 
          [&](double x) {return this->radial_contr.compute_neighbour_contribution(x, 0.5);};
      intp.initialize(func, x1, x2, error_bound); 
      
      int matrix_size = max_radial*(max_angular+1);
      Matrix_t intp_val = Matrix_t::Zero(ref_points.size(), matrix_size);
      Matrix_t intp_ref = Matrix_t::Zero(ref_points.size(), matrix_size);
      for (int i{0}; i<ref_points.size()-1; i++) {
        intp_val.row(i) = Eigen::Map<Vector_t>(intp.interpolate(ref_points(i)).data(),matrix_size);
        intp_ref.row(i) = Eigen::Map<Vector_t>(func(ref_points(i)).data(),matrix_size);
      }
      double error{(intp_val-intp_ref).array().abs().colwise().mean().maxCoeff()};
      BOOST_CHECK_LE(error, error_bound);
    }

  BOOST_AUTO_TEST_SUITE_END();
  }  // namespace math 
}  // namespace rascal
