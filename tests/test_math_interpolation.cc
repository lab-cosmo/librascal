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

  BOOST_AUTO_TEST_SUITE(MathInterpolatorTests);

  using HuntInterpolator = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Hunt>
        >;
  using UniformInterpolator = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;
  using interpolator_fixtures = boost::mpl::list<InterpolatorFixture<HuntInterpolator>,
                   InterpolatorFixture<UniformInterpolator>>;                     

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(math_interpolator_tests, Fix,
                                   interpolator_fixtures, Fix) {
    auto intp{Fix::intp};
    auto identity_func{Fix::identity_func};
    auto exp_func{Fix::exp_func};
    auto mean_error_bound{Fix::mean_error_bound};

    Vector_t ref_points = Vector_t::LinSpaced(Fix::nb_ref_points, Fix::x1, Fix::x2); 
    double error, intp_val, intp_ref;
    std::function<double(double)> func;

    // TODO(alex) make function array
    // Test for identity function
    func = identity_func;
    intp.initialize(func, Fix::x1, Fix::x2, Fix::mean_error_bound);
    for (int i{0}; i<ref_points.size()-1; i++) {
      intp_val = intp.interpolate(ref_points(i));
      intp_ref = func(ref_points(i));
      error = std::abs(intp_val - intp_ref);
      BOOST_CHECK_LE(error, mean_error_bound);
    }

    // Test for exp function
    func = exp_func;
    intp.initialize(func, Fix::x1, Fix::x2, Fix::mean_error_bound);
    for (int i{0}; i<ref_points.size()-1; i++) {
      intp_val = intp.interpolate(ref_points(i));
      intp_ref = func(ref_points(i));
      error = std::abs(intp_val - intp_ref);
      BOOST_CHECK_LE(error, mean_error_bound);
    }
    
    func  = [&](double x) {return this->hyp1f1.calc(x);};
    //func = std::bind(hyp1f1_function_generator, a, b, std::placeholders::_1);
    intp.initialize(func, Fix::x1, Fix::x2, Fix::mean_error_bound);
    for (int i{0}; i<ref_points.size()-1; i++) {
      intp_val = intp.interpolate(ref_points(i));
      intp_ref = func(ref_points(i));
      error = std::abs(intp_val - intp_ref);
      BOOST_CHECK_LE(error, mean_error_bound);
    }
  }
    
  // Test for identity function
  BOOST_FIXTURE_TEST_CASE(math_interpolator_identity_test, InterpolatorFixture<UniformInterpolator>) {
    Vector_t ref_points = Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    this->intp.initialize(this->identity_func, this->x1, this->x2, this->mean_error_bound); 
    BOOST_CHECK_LE(intp_ref_mean_error(this->intp, ref_points), this->mean_error_bound);
  }

  // Test for identity function
  BOOST_FIXTURE_TEST_CASE(math_interpolator_exp_test, InterpolatorFixture<UniformInterpolator>) {
    Vector_t ref_points = Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    this->intp.initialize(this->exp_func, this->x1, this->x2, this->mean_error_bound); 
    BOOST_CHECK_LE(intp_ref_mean_error(this->intp, ref_points), this->mean_error_bound);
  }
  // tests if functions are well approximated up to precisicion
  BOOST_FIXTURE_TEST_CASE(math_interpolator_hyp1f1_test, InterpolatorFixture<UniformInterpolator>) {
    Vector_t ref_points = Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    std::function<double(double)> func = [&](double x) {return this->hyp1f1.calc(x);};
    //std::function<double(double)> func = std::bind(hyp1f1_function_generator, a, b, std::placeholders::_1);
    intp.initialize(func, this->x1, this->x2, this->mean_error_bound); 
    BOOST_CHECK_LE(intp_ref_mean_error(this->intp, ref_points), this->mean_error_bound);
    //  }
    //}
  }

  BOOST_FIXTURE_TEST_CASE(math_interpolator_radial_contribution_test, InterpolatorFixture<UniformInterpolator>) {
    Vector_t ref_points = Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
    std::function<double(double)> func =
        [&](double x) {return this->radial_contr.compute_contribution<rascal::internal::AtomicSmearingType::Constant>(x, 0.5)(0,0);};
    intp.initialize(func, this->x1, this->x2, this->mean_error_bound); 
    BOOST_CHECK_LE(intp_ref_mean_error(this->intp, ref_points), this->mean_error_bound);
  }

  using DefaultInterpolatorVectorized = InterpolatorVectorized<
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;

  BOOST_FIXTURE_TEST_CASE(math_interpolator_vectorized_test, InterpolatorFixture<DefaultInterpolatorVectorized>) {
      Vector_t ref_points = Vector_t::LinSpaced(nb_ref_points, x1, x2); 
      double error;

      // test scalar interpolator and vectorized interpolator give same results

      auto intp_scalar = UniformInterpolator();
      std::function<double(double)> func_scalar =
          [&](double x) {return this->radial_contr.compute_contribution<rascal::internal::AtomicSmearingType::Constant>(x, 0.5)(0,0);};
      intp_scalar.initialize(func_scalar, x1, x2, mean_error_bound);

      std::function<Matrix_t(double)> func = [&](double x) {return this->radial_contr.compute_contribution<rascal::internal::AtomicSmearingType::Constant>(x, 0.5);};
      ////// change to function generator
      ////func = std::bind(hyp1f1_function_generator, a, b, std::placeholders::_1);
      //intp.initialize(func, 0,1, mean_error_bound); 
      intp.initialize(func, x1, x2, mean_error_bound, 512); 


      // Checks if vectorized interpolator and scalar interpolator return same results 
      double intp_vec_val, intp_scalar_val;
      for (int i{0}; i<ref_points.size()-1; i++) {
        intp_vec_val = intp.interpolate(ref_points(i))(0,0);
        intp_scalar_val = intp_scalar.interpolate(ref_points(i)); 
        error = std::abs(intp_vec_val - intp_scalar_val);
        BOOST_CHECK_LE(error, 1e-10);// TODO(alex) use global tolerance
      }


      // checks if error bound is fulfilled
      //Matrix_t intp_val = Matrix_t::Zero(max_radial,max_angular);
      //Matrix_t intp_ref = Matrix_t::Zero(max_radial,max_angular);
      // Test for max error function
      // TODO(alex) find a better way of testing this, error function is changed
      //for (int i{0}; i<ref_points.size()-1; i++) {
      //  intp_val = intp.interpolate(ref_points(i));
      //  intp_ref = func(ref_points(i));
      //  error = (intp_val - intp_ref).array().abs().maxCoeff();
      //  BOOST_CHECK_LE(error, 1e-6);
      //}
      error = 0.0;
      int matrix_size = max_radial*max_angular;
      Matrix_t intp_val = Matrix_t::Zero(ref_points.size(), matrix_size);
      Matrix_t intp_ref = Matrix_t::Zero(ref_points.size(), matrix_size);
      for (int i{0}; i<ref_points.size()-1; i++) {
        intp_val.row(i) = Eigen::Map<Vector_t>(func(ref_points(i)).data(),matrix_size);
        intp_ref.row(i) = Eigen::Map<Vector_t>(intp.interpolate(ref_points(i)).data(),matrix_size);
      }
      error = (intp_val-intp_ref).array().abs().colwise().mean().maxCoeff();
      BOOST_CHECK_LE(error, mean_error_bound);
    }

  BOOST_AUTO_TEST_SUITE_END();
  }  // namespace math 
}  // namespace rascal
