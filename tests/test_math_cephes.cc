/**
 * file   test_math_cephes.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief Test interface between librascal and cephes
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_math.hh"

namespace rascal {
  constexpr static double math_tol{1e-14};

  BOOST_AUTO_TEST_SUITE(MathCEPHESTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(math_hyp2f1_test, ManagerFixture_math) {
    for (int ii{0}; ii < 3; ++ii) {
      double val{math::hyp2f1(numbers(0, ii), numbers(1, ii),
                              numbers(2, ii), numbers(3, ii))};
      if (vebose) {
        std::cout << std::setprecision(14) << val <<", ";
      }
      auto error{std::abs(val-results_hyp2f1[ii])};
      BOOST_CHECK_LE(error, math_tol);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
