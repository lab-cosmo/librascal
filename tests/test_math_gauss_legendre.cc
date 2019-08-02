/**
 * file   test_math_gauss_legendre.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   21 May 2019
 *
 * @brief Test the implementation of gauss legendre quadrature against numpy
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

  BOOST_AUTO_TEST_SUITE(MathGaussLegendreTests);

  /* ---------------------------------------------------------------------- */
  /**
   * Check the implementation of hyp1f1 against mpmath v1.1.0
   */
  BOOST_FIXTURE_TEST_CASE(math_hyp1f1_test, GaussLegendreRefFixture) {
    for (auto & data : this->ref_data) {
      double a{data["a"]}, b{data["b"]};
      int order{data["order"]};
      auto points_ref{data["points"].get<std::vector<double>>()};
      auto weights_ref{data["weights"].get<std::vector<double>>()};

      auto val{math::compute_gauss_legendre_points_weights(a, b, order)};

      for (int ii{0}; ii < order; ++ii) {
        double points_rel_error{(val(ii, 0) - points_ref[ii])};
        double weigths_rel_error{(val(ii, 1) - weights_ref[ii])};

        if ((points_rel_error > math::dbl_ftol or
             weigths_rel_error > math::dbl_ftol) and
            this->verbose) {
          std::cout << " a=" << a << " b=" << b << " order=" << order
                    << " point_err=" << points_rel_error
                    << " weight_err=" << weigths_rel_error << std::endl;
        }

        BOOST_CHECK_LE(points_rel_error, math::dbl_ftol);
        BOOST_CHECK_LE(weigths_rel_error, math::dbl_ftol);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
