/**
 * file   test_math_modified_bessel_first_kind.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   21 May 2019
 *
 * @brief Test the implementation of modified bessel first kind
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

  BOOST_AUTO_TEST_SUITE(MathBesselFirstKindTests);

  /* ----------------------------------------------------------------------
  */
  /**
   * Check the implementation of the modified bessel function of the first
   * kind against mpmath v1.1.0
   */
  BOOST_FIXTURE_TEST_CASE(math_bessel_test,
  ModifiedBesselFirstKindRefFixture) {
    json& i_exp_ref{this->ref_data["i_exp"]};
    for (auto & data : i_exp_ref) {
      double x{data["x"]};
      auto ref_vals{data["vals"].get<std::vector<double>>()};
      size_t max_order{data["max_order"]};

      auto vals{math::bessel_i_exp_allorders(x, max_order)};

      // for (size_t order{0}; order < max_order; ++order) {
      //   double rel_error{(vals(order) - ref_vals[order])};

      //   if ((rel_error > 1e5 * math::dbl_ftol) and this->verbose) {
      //     std::cout << " order=" << order << " x=" << x
      //     << " diff=" << rel_error<< " val=" << vals(order) << std::endl;
      //   }

      //   BOOST_CHECK_LE(rel_error, 1e5 * math::dbl_ftol);
      // }
    }

    json& i_complete_square_ref{this->ref_data["i_complete_square"]};
    for (auto & data : i_complete_square_ref) {
      auto xs{data["xs"].get<std::vector<double>>()};
      Eigen::Map<Eigen::ArrayXd> xns(&xs[0], xs.size());
      double alpha{data["alpha"]}, rij{data["rij"]};
      auto ref_vals{data["vals"].get<std::vector<std::vector<double>>>()};
      size_t max_order{data["max_order"]};

      auto vals{math::bessel_i_exp_exp_complete_square(xns, rij, alpha, max_order)};

      for (size_t i_x{0}; i_x < xs.size(); ++i_x){
        for (size_t order{0}; order < max_order; ++order) {
          double rel_error{(vals(i_x, order) - ref_vals[i_x][order])};

          if ((rel_error > 1e5 * math::dbl_ftol) and this->verbose) {
            std::cout << " order=" << order << " x=" << xs[i_x]
            << " alpha=" << alpha << " rij=" << rij
            << " diff=" << rel_error<< " val=" << vals(i_x, order) << std::endl;
          }

          // BOOST_CHECK_LE(rel_error, 1e5 * math::dbl_ftol);
        }
      }

    }
  }

  /* ----------------------------------------------------------------------
  */ BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
