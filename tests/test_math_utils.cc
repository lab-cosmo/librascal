/**
 * file   test_math_math.cc
 *
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   23 November 2018
 *
 * @brief Test own math utilities (spherical harmonics)
 *
 * Copyright © 2018  Max Veit, COSMO (EPFL), LAMMM (EPFL)
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
  constexpr static double math_tol{1e-14};

  BOOST_AUTO_TEST_SUITE(MathUtilsTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test,
                          SphericalHarmonicsRefFixture) {
    for (size_t vec_idx{0}; vec_idx < unit_vectors.size(); vec_idx++) {
      Eigen::Vector3d direction(unit_vectors[vec_idx].data());
      size_t max_angular = harmonics[vec_idx].size() - 1;
      Eigen::MatrixXd sph_harm(max_angular+1, 2*max_angular+1);
      sph_harm = math::compute_spherical_harmonics(
          direction, max_angular, sph_harm);
      for (size_t ang_idx{0}; ang_idx < max_angular + 1; ang_idx++) {
        for (size_t m_idx{0}; m_idx < ang_idx + 2; m_idx++) {
          // Check harmonics[vec_idx][ang_idx][m_idx] against
          // sph_harm[ang_idx][m_idx]
        }
      }
    }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
