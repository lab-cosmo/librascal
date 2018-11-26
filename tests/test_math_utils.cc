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
  // Double ftol defined in test_math.hh (currently 100*eps, so about the same
  // as below)
  //constexpr static double math_tol{1e-14};

  BOOST_AUTO_TEST_SUITE(MathUtilsTests);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test,
                          SphericalHarmonicsRefFixture) {
    for (size_t vec_idx{0}; vec_idx < unit_vectors.size(); vec_idx++) {
      Eigen::Vector3d direction(unit_vectors[vec_idx].data());
      size_t max_angular = harmonics[vec_idx].size() - 1;
      Eigen::MatrixXd computed_harmonics = math::compute_spherical_harmonics(
          direction, max_angular);
      if (verbose) {
        std::cout << "Testing unit vector: " << direction << std::endl;
        std::cout << "Max angular momentum: l_max=" << max_angular << std::endl;
        std::cout << "Computed harmonics size: " << computed_harmonics.rows();
        std::cout << " by " << computed_harmonics.cols() << std::endl;
      }
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        if (verbose) {
          std::cout << std::setprecision(10) << "Coefficients for l=";
          std::cout << angular_l << ": ";
          std::cout << computed_harmonics.block(
                          angular_l, 0, 1, 2*angular_l + 1);
          std::cout << std::endl;
        }
        for (size_t m_idx{0}; m_idx < 2*angular_l + 1; m_idx++) {
          // Check both the harmonics and their order in memory
          auto error{std::abs(computed_harmonics(angular_l, m_idx)
                              - harmonics[vec_idx][angular_l][m_idx])};

          BOOST_CHECK_LE(error, math::dbl_ftol);
        }
      }
      // TODO(max-veit) do we make sure that the rest of the array is zero?
    }
  }

  BOOST_FIXTURE_TEST_CASE(math_alp_test, SphericalHarmonicsRefFixture) {
    for (size_t vec_idx{0}; vec_idx < unit_vectors.size(); vec_idx++) {
      Eigen::Vector3d direction(unit_vectors[vec_idx].data());
      size_t max_angular = harmonics[vec_idx].size() - 1;
      Eigen::MatrixXd computed_alps = math::compute_assoc_legendre_polynom(
          direction[2], max_angular);
      if (verbose) {
        std::cout << "Testing unit vector: " << direction.transpose() << std::endl;
        std::cout << "Max angular momentum: l_max=" << max_angular << std::endl;
        std::cout << "Computed harmonics size: " << computed_alps.rows();
        std::cout << " by " << computed_alps.cols() << std::endl;
      }
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        if (verbose) {
          std::cout << std::setprecision(10) << "Computed ALPs for l=";
          std::cout << angular_l << ": ";
          std::cout << computed_alps.block(angular_l, 0, 1, angular_l + 1);
          std::cout << std::endl;
        }
        for (size_t m_idx{0}; m_idx < angular_l + 1; m_idx++) {
          // Check both the ALPs and their order in memory
          auto error{std::abs(computed_alps(angular_l, m_idx)
                              - alps[vec_idx][angular_l][m_idx])};
          BOOST_CHECK_LE(error, math::dbl_ftol);
        }
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE(math_cos_sin_mphi_test,
                          SphericalHarmonicsRefFixture) {
    size_t max_m = 10;
    Eigen::VectorXd phi_test(5);
    phi_test << 0, 0.1, math::PI/4, math::PI, math::PI*2;
    for (size_t phi_idx{0}; phi_idx < 5; phi_idx++) {
      double cos_phi = std::cos(phi_test(phi_idx));
      double sin_phi = std::sin(phi_test(phi_idx));
      Eigen::MatrixXd cos_sin_m_phi = math::compute_cos_sin_angle_multiples(
          cos_phi, sin_phi, max_m);
      if(verbose) {
        std::cout << "Cos | sin (mφ), φ=" << phi_test(phi_idx) << std::endl;
        std::cout << cos_sin_m_phi.transpose() << std::endl;
      }
      for (size_t m_idx{0}; m_idx < max_m; m_idx++) {
        auto cos_error{std::abs(cos_sin_m_phi(m_idx, 0) -
                                std::cos(m_idx*phi_test(phi_idx)))};
        auto sin_error{std::abs(cos_sin_m_phi(m_idx, 1) -
                                std::sin(m_idx*phi_test(phi_idx)))};
        BOOST_CHECK_LE(cos_error, math::dbl_ftol);
        BOOST_CHECK_LE(sin_error, math::dbl_ftol);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
