/**
 * file   test_math_math.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   21 May 2019
 *
 * @brief Test the implementation of Hyp1f1 against mpmath
 *
 * Copyright © 2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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


  BOOST_AUTO_TEST_SUITE(MathHyp1f1Tests);

  /* ---------------------------------------------------------------------- */
  /**
   * Check the implementation of hyp1f1 against mpmath v1.1.0
   */
  BOOST_FIXTURE_TEST_CASE(math_hyp1f1_test,
                          Hyp1F1RefFixture) {
    for (auto& data : this->ref_data) {
      double a{data["a"]},b{data["b"]},z{data["z"]},hyp1f1_ref{data["val"]},hyp1f1_der_ref{data["der"]};
      math::Hyp1f1 func{a,b,200,1e-13};
      double val{func.calc(z)};
      double der{func.calc(z, true)};
      // TODO(felix) find a proper scheme for numerical derivatives
      // double h{1e-9};
      // // centered finite difference
      // // double hyp1f1_num_der{(func.calc(z+h)-func.calc(z-h)) / (2*h)};
      // double hyp1f1_num_der{
      //   (-func.calc(z+2*h)+8*func.calc(z+h)-8*func.calc(z-h)+func.calc(z-2*h)) / (12*h)};

      double rel_error{std::abs((hyp1f1_ref-val)/hyp1f1_ref)};
      if (rel_error > 10*math::dbl_ftol and this->verbose) {
        std::cout << " a=" << a<< " b=" << b<< " z=" << z<< " ref=" << hyp1f1_ref<< " impl="<< val<< " z_switch=";
        std::cout << func.z_asympt << std::endl;
      }
      BOOST_CHECK_LE(rel_error, 10*math::dbl_ftol);

      double rel_der_error{std::abs((hyp1f1_der_ref-der)/hyp1f1_der_ref)};
      if (rel_der_error > 10*math::dbl_ftol and this->verbose) {
        std::cout << "Derivative a=" << a<< " b=" << b<< " z=" << z<< " ref=" << hyp1f1_der_ref<< " impl="<< der<< " z_switch=";
        std::cout << func.z_asympt << std::endl;
      }
      BOOST_CHECK_LE(rel_der_error, 10*math::dbl_ftol);

      // double der_consistency_rel_error{std::abs((hyp1f1_num_der-der)/hyp1f1_num_der)};
      // if (der_consistency_rel_error > 1000*math::dbl_ftol and this->verbose) {
      //   std::cout << "Derivative consistency a=" << a<< " b=" << b<< " z=" << z<< " num_der=" << hyp1f1_num_der<< " impl="<< der << " rel_diff="<<der_consistency_rel_error<< " z_switch=";
      //   std::cout << func.z_asympt << std::endl;
      // }
      // // BOOST_CHECK_LE(der_consistency_rel_error, 10*math::dbl_ftol);
    }
  }


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
