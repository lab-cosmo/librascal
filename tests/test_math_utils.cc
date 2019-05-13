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
      Eigen::MatrixXd computed_alps =
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

  constexpr size_t test_max_angular = 3;
  using function_list = boost::mpl::list<
    SphericalHarmonicsWithGradients<test_max_angular>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(math_gradients_test, Function_provider_t,
                                   function_list, GradientTestFixture) {
    Function_provider_t function_calculator{};
    Eigen::MatrixXd values;
    Eigen::MatrixXd jacobian;
    Eigen::RowVectorXd direction_vector;
    Eigen::VectorXd displacement_direction;
    Eigen::Vector3d displacement;
    Eigen::MatrixXd directional;
    Eigen::MatrixXd fd_derivatives;
    Eigen::MatrixXd fd_error_cwise;
    // This error isn't going to be arbitrarily small, due to the interaction of
    // finite-difference and finite precision effects.  Just set it to something
    // reasonable and check it explicitly if you really want to be sure (paying
    // attention to the change of the finite-difference gradients from one step
    // to the next).  The automated test is really more intended to be a sanity
    // check on the implementation anyway.
    constexpr double fd_error_tol = 1E-8;
    for (auto inputs_it{function_inputs.begin()};
         inputs_it != function_inputs.end(); inputs_it++) {
      direction_vector = Eigen::Map<Eigen::RowVectorXd>
                                          (inputs_it->data(), 1, n_arguments);
      values = function_calculator.f(direction_vector);
      jacobian = function_calculator.grad_f(direction_vector);
      if (verbose) {
        std::cout << "Direction vector: " << direction_vector << std::endl;
        std::cout << "Values:" << values << std::endl;
        std::cout << "Jacobian:" << jacobian << std::endl;
      }
      for (int disp_idx{0}; disp_idx < displacement_directions.rows();
           disp_idx++) {
        displacement_direction = displacement_directions.row(disp_idx);
        // Compute the directional derivative(s)
        directional = displacement_direction.adjoint() * jacobian;
        if (verbose) {
          std::cout << "FD direction: " << displacement_direction << std::endl;
          std::cout << "Analytical derivative: " << directional << std::endl;
        }
        double min_error{HUGE_VAL};
        Eigen::MatrixXd fd_last{Eigen::MatrixXd::Zero(1, directional.size())};
        for (double dx = 1E-2; dx > 1E-10; dx *= 0.1) {
          std::cout << "dx = " << dx << "\t";
          displacement = dx * displacement_direction;
          // Compute the finite-difference derivative using a
          // centred-difference approach
          fd_derivatives = 0.5 / dx * (
            function_calculator.f(direction_vector + displacement.adjoint())
          - function_calculator.f(direction_vector - displacement.adjoint()));
          double fd_error{0.};
          double fd_quotient{0.};
          size_t nonzero_count{0};
          for (int dim_idx{0}; dim_idx < fd_derivatives.size(); dim_idx++) {
            if (std::abs(directional(dim_idx)) < 10*math::dbl_ftol) {
              fd_error += fd_derivatives(dim_idx);
            } else {
              fd_quotient += (fd_derivatives(dim_idx) / directional(dim_idx));
              fd_error += (fd_derivatives(dim_idx) - directional(dim_idx))
                          / directional(dim_idx);
              ++nonzero_count;
            }
          }
          if (nonzero_count > 0) {
            fd_quotient = fd_quotient / nonzero_count;
          }
          fd_error = fd_error / fd_derivatives.size();
          std::cout << "Average rel FD error: " << fd_error << "\t";
          std::cout << "Average FD quotient:  " << fd_quotient << std::endl;
          if (std::abs(fd_error) < min_error) {min_error = std::abs(fd_error);}
          if (verbose) {
            fd_error_cwise = (fd_derivatives - directional);
            std::cout << "error            = " << fd_error_cwise << std::endl;
            std::cout << "(FD derivative   = " << fd_derivatives << ")";
            std::cout << std::endl;
            std::cout << "(minus last step = " << fd_derivatives - fd_last;
            std::cout << ")" << std::endl;
          }
          fd_last = fd_derivatives;
        } // for (double dx...) (displacement magnitudes)
        BOOST_CHECK_SMALL(min_error, fd_error_tol);
      } // for (int disp_idx...) (displacement directions)
    } // for (auto inputs_it...) (function inputs)
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
