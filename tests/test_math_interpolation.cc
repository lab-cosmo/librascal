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

  // tests if functions are well approximated up to precisicion
  BOOST_FIXTURE_TEST_CASE(math_interpolator_test, InterpolatorFixture) {
    auto intp{Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::HeapBased>,
      SearchMethod<SearchMethod_t::Hunt>
        >()};
    Vector_t points(11);
    points << 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.;
    double error, intp_val, intp_ref;
    std::function<double(double)> func;

    // TODO(alex) make function array
    // Test for identity function
    func = this->identity_func;
    intp.initalize(func, 0,1, this->precision); 
    for (int i{0}; i<points.size()-1; i++) {
      intp_val = intp.interpolate(points(i));
      intp_ref = func(points(i));
      error = std::abs(intp_val - intp_ref);
      BOOST_CHECK_LE(error, this->precision);
    }

    // Test for exp function
    func = this->exp_func;
    intp.initalize(func, 0,1, this->precision); 
    for (int i{0}; i<points.size()-1; i++) {
      intp_val = intp.interpolate(points(i));
      intp_ref = func(points(i));
      error = std::abs(intp_val - intp_ref);
      BOOST_CHECK_LE(error, this->precision);
    }
    
    size_t max_angular = 20;
    size_t max_radial = 20;
    double a,b;

    bool verbose{false};
    size_t l{max_angular};
    size_t n{max_radial};
    // TODO(alex) this test would take 5 seconds, therefore
    // I only test subset
    //for (size_t l{0}; l<max_angular; l++) {
    //  for (size_t n{0}; n<max_radial; n++) {
        if (verbose) {
          std::cout << "Testing for n="<<n<<" and l="<<l<<std::endl;
        }
        a = 0.5*(n+l+3);
        b = l+1.5;
        math::Hyp1f1 hyp1f1{a, b, 200, 1e-15};
        func  = [&hyp1f1](double x) {return hyp1f1.calc(x);};
        //func = std::bind(hyp1f1_function_generator, a, b, std::placeholders::_1);
        intp.initalize(func, 0,1, this->precision); 
        for (int i{0}; i<points.size()-1; i++) {
          intp_val = intp.interpolate(points(i));
          intp_ref = func(points(i));
          error = std::abs(intp_val - intp_ref);
          BOOST_CHECK_LE(error, this->precision);
        }
    //  }
    //}
  }

  BOOST_AUTO_TEST_SUITE_END();
  }  // namespace math 
}  // namespace rascal
