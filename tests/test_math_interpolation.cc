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
  }

  BOOST_AUTO_TEST_SUITE_END();
  }  // namespace math 
}  // namespace rascal