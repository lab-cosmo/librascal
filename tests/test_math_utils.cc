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
  // constexpr static double math_tol{1e-14};

  BOOST_AUTO_TEST_SUITE(MathUtilsTests);

  /* ---------------------------------------------------------------------- */
  /**
   * Check the spherical harmonics implementation
   *
   * The spherical harmonics are checked against the equivalents computed using
   * SciPy, adjusted to use our normalization for the real harmonics.  They can
   * be generated with the following Python code:
   *
   * ```python
   * import numpy as np
   * from scipy import special as spf
   * harmonics = np.zeros((lmax+1, 2*lmax + 1, unit_vectors.shape[0]))
   * for l in range(lmax+1):
   *     harmonics[l, l, :] = spf.sph_harm(0, l, phis, thetas)
   *     for m in range(1, l+1):
   *         complex_harmonics = spf.sph_harm(m, l, phis, thetas)
   *         harmonics[l, l+m, :] = np.real(complex_harmonics)*np.sqrt(2)
   *         harmonics[l, l-m, :] = -1*np.imag(complex_harmonics)*np.sqrt(2)
   * ```
   *
   * where `lmax` denotes the maximum angular momentum number to compute them
   * up to, and `thetas` and `phis` are lists (arrays) of angles at which to
   * compute the harmonics; `theta` is the angle off the z-axis and `phi` is the
   * angle (projected into the x-y plane) off the x-axis.
   *
   * The code above works with SciPy v1.1.0.
   */
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test,
                          SphericalHarmonicsRefFixture) {
    for (size_t vec_idx{0}; vec_idx < unit_vectors.size(); vec_idx++) {
      Eigen::Vector3d direction(unit_vectors[vec_idx].data());
      size_t max_angular = harmonics[vec_idx].size() - 1;
      Eigen::VectorXd computed_harmonics =
          math::compute_spherical_harmonics(direction, max_angular);
      if (verbose) {
        std::cout << "Testing unit vector: " << direction << std::endl;
        std::cout << "Max angular momentum: l_max=" << max_angular << std::endl;
        std::cout << "Number of computed harmonics: ";
        std::cout << computed_harmonics.size() << std::endl;
      }
      size_t lm_collective_idx{0};
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        if (verbose) {
          std::cout << std::setprecision(10) << "Coefficients for l=";
          std::cout << angular_l << ": ";
          std::cout << computed_harmonics.segment(lm_collective_idx,
                                                  2*angular_l + 1);
          std::cout << std::endl;
        }
        for (size_t m_idx{0}; m_idx < 2 * angular_l + 1; m_idx++) {
          // Check both the harmonics and their order in memory
          auto error{std::abs(computed_harmonics(lm_collective_idx + m_idx) -
                              harmonics[vec_idx][angular_l][m_idx])};
          BOOST_CHECK_LE(error, math::dbl_ftol);
        }
        lm_collective_idx += (2*angular_l + 1);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Check the Associated Legendre Polynomials recurrent implementation
   *
   * The ALPs are checked against the equivalents computed using
   * SciPy, adjusted to use our normalization (see
   * https://arxiv.org/pdf/1410.1748.pdf; these functions are $\bar{P}_l^m$ in
   * their notation).  They can be generated using the following Python code:
   *
   * ```python
   * import numpy as np
   * from scipy import special as spf
   * alp_normfacts = np.zeros((lmax+1, lmax+1))
   * for l in range(lmax+1):
   *     for m in range(l+1):
   *         alp_normfacts[l, m] = np.sqrt(
   *             (2*l + 1)/(2*np.pi) / np.prod(np.arange(l-m+1, l+m+1)))
   * for d_idx in range(n_directions):
   *     alps[:, :, d_idx] = (spf.lpmn(l, l, np.cos(thetas[d_idx]))[0].T
   *                          * alp_normfacts[:, :, np.newaxis])
   * ```
   *
   * where `lmax` and `thetas` are defined as above and the matrix product at
   * the end is taken elementwise.
   */
  BOOST_FIXTURE_TEST_CASE(math_associated_legendre_polynomial_test,
                          SphericalHarmonicsRefFixture) {
    for (size_t vec_idx{0}; vec_idx < unit_vectors.size(); vec_idx++) {
      Eigen::Vector3d direction(unit_vectors[vec_idx].data());
      size_t max_angular = harmonics[vec_idx].size() - 1;
      math::Array_t computed_alps =
          math::compute_assoc_legendre_polynom(direction[2], max_angular);
      if (verbose) {
        std::cout << "Testing unit vector: " << direction.transpose();
        std::cout << std::endl;
        std::cout << "Max angular momentum: l_max=" << max_angular << std::endl;
        std::cout << "Computed harmonics size: " << computed_alps.rows();
        std::cout << " by " << computed_alps.cols() << std::endl;
      }
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        if (verbose) {
          std::cout << std::setprecision(10) << "Computed ALPs for l=";
          std::cout << angular_l << ": ";
          std::cout << computed_alps.block(angular_l, 0, 1, angular_l + 2);
          std::cout << std::endl;
        }
        for (size_t m_idx{0}; m_idx < angular_l + 1; m_idx++) {
          // Check both the ALPs and their order in memory
          auto error{std::abs(computed_alps(angular_l, m_idx) -
                              alps[vec_idx][angular_l][m_idx])};
          BOOST_CHECK_LE(error, math::dbl_ftol);
        }
        for (size_t m_idx{angular_l + 1}; m_idx < max_angular + 2; m_idx++) {
          BOOST_CHECK_EQUAL(computed_alps(angular_l, m_idx), 0.);
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(math_cos_sin_mphi_test,
                          SphericalHarmonicsRefFixture) {
    size_t max_m = 10;
    Eigen::VectorXd phi_test(5);
    phi_test << 0, 0.1, math::PI / 4, math::PI, math::PI * 2;
    for (size_t phi_idx{0}; phi_idx < 5; phi_idx++) {
      double cos_phi = std::cos(phi_test(phi_idx));
      double sin_phi = std::sin(phi_test(phi_idx));
      Eigen::MatrixXd cos_sin_m_phi =
          math::compute_cos_sin_angle_multiples(cos_phi, sin_phi, max_m);
      if (verbose) {
        std::cout << "Cos | sin (mφ), φ=" << phi_test(phi_idx) << std::endl;
        std::cout << cos_sin_m_phi.transpose() << std::endl;
      }
      for (size_t m_idx{0}; m_idx < max_m; m_idx++) {
        auto cos_error{std::abs(cos_sin_m_phi(m_idx, 0) -
                                std::cos(m_idx * phi_test(phi_idx)))};
        auto sin_error{std::abs(cos_sin_m_phi(m_idx, 1) -
                                std::sin(m_idx * phi_test(phi_idx)))};
        BOOST_CHECK_LE(cos_error, math::dbl_ftol);
        BOOST_CHECK_LE(sin_error, math::dbl_ftol);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(spherical_harmonics_gradient_test) {
    constexpr size_t test_max_angular = 3;
    test_gradients<SphericalHarmonicsWithGradients<test_max_angular>>(
      "reference_data/spherical_harmonics_gradient_test.json");
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
