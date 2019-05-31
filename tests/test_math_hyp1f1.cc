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
  BOOST_FIXTURE_TEST_CASE(math_hyp1f1_test, Hyp1F1RefFixture) {
    for (auto & data : this->ref_data) {
      double a{data["a"]}, b{data["b"]}, z{data["z"]}, hyp1f1_ref{data["val"]},
          hyp1f1_der_ref{data["der"]};
      math::Hyp1f1 func{a, b, 200, 1e-13};
      double val{func.calc(z)};
      double der{func.calc(z, true)};

      double h{1e-5};
      double hyp1f1_num_der{func.calc_numerical_derivative(z, h)};

      // check if hyp1f1 is consistent with the
      // mpmath reference
      double rel_error{std::abs((hyp1f1_ref - val) / hyp1f1_ref)};
      if (rel_error > 15 * math::dbl_ftol and this->verbose) {
        std::cout << " a=" << a << " b=" << b << " z=" << z
                  << " ref=" << hyp1f1_ref << " impl=" << val << " z_switch=";
        std::cout << func.z_asympt << std::endl;
      }
      BOOST_CHECK_LE(rel_error, 15 * math::dbl_ftol);

      // check if the analytical derivatives are consistent with the
      // mpmath reference
      double rel_der_error{std::abs((hyp1f1_der_ref - der) / hyp1f1_der_ref)};
      if (rel_der_error > 15 * math::dbl_ftol and this->verbose) {
        std::cout << "Derivative a=" << a << " b=" << b << " z=" << z
                  << " ref=" << hyp1f1_der_ref << " impl=" << der
                  << " rel_err=" << rel_der_error << " z_switch=";
        std::cout << func.z_asympt << std::endl;
      }
      BOOST_CHECK_LE(rel_der_error, 15 * math::dbl_ftol);

      // check if the numerical derivatives are consistent with the
      // analytical ones
      double der_consistency_rel_error{
          std::abs((hyp1f1_num_der - der) / hyp1f1_num_der)};
      if (der_consistency_rel_error > 1e5 * math::dbl_ftol and this->verbose) {
        std::cout << "Derivative consistency a=" << a << " b=" << b
                  << " z=" << z << " num_der=" << hyp1f1_num_der
                  << " impl=" << der
                  << " rel_diff=" << der_consistency_rel_error << " z_switch=";
        std::cout << func.z_asympt << std::endl;
      }

      BOOST_CHECK_LE(der_consistency_rel_error, 1e5 * math::dbl_ftol);
    }
  }

  BOOST_FIXTURE_TEST_CASE(math_hyp1f1_spherical_expansion_test,
                          Hyp1f1SphericalExpansionFixture) {
    for (size_t i_rc{0}; i_rc < this->rcs.size(); ++i_rc) {
      auto & rc{this->rcs[i_rc]};
      auto & fac_b{this->facs_b[i_rc]};
      for (auto & r_ij : this->r_ijs) {
        if (r_ij >= rc) {
          continue;
        }
        for (auto & fac_a : this->fac_as) {
          for (size_t ii{0}; ii < this->hyp1f1.size(); ++ii) {
            hyp1f1[ii].calc(r_ij, fac_a, fac_b[ii], true);
            hyp1f1_recursion[ii].calc(r_ij, fac_a, fac_b[ii], true);
            auto hyp1f1_val{hyp1f1[ii].get_values()};
            auto hyp1f1_recursion_val{hyp1f1_recursion[ii].get_values()};
            auto diff_val{(hyp1f1_val - hyp1f1_recursion_val).array().abs() /
                          hyp1f1_val.array()};
            auto hyp1f1_der{hyp1f1[ii].get_derivatives()};
            auto hyp1f1_recursion_der{hyp1f1_recursion[ii].get_derivatives()};
            auto diff_der{(hyp1f1_der - hyp1f1_recursion_der).array().abs() /
                          hyp1f1_der.array()};
            BOOST_CHECK_LE(diff_val.mean(), 3 * math::dbl_ftol);
            BOOST_CHECK_LE(diff_der.mean(), 3 * math::dbl_ftol);
            if (verbose) {
              std::cout << "diff_val= " << diff_val.mean()
                        << " diff_der=" << diff_der.mean() << std::endl;
              for (int ii{0}; ii < diff_val.rows(); ++ii) {
                for (int jj{0}; jj < diff_val.cols(); ++jj) {
                  std::cout << diff_val(ii, jj) << ", ";
                }
              }
            }
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
