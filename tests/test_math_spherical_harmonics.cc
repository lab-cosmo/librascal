/**
 * file   test_math_spherical_harmonics.cc
 *
 * @section LICENSE
 *
 * Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_math.hh"

#include <boost/test/unit_test.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(MathSphericalHarmonicsTests);

  // We assume that max_angular_l does not change in the reference data.
  // TODO(alex) make this more clear by changing the data structure of the
  // reference data
  // TODO(alex) replace info and verbose usage of VerbosityValue NORMAL DEBUG
  /* Rescued from `test_math_utils.cc`: Keep or refer to the script used to
   * generate reference values
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
   */
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test,
                          SphericalHarmonicsClassRefFixture) {
    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test " << std::endl;
      std::cout << ">> with number_of_unit_vectors=" << this->ref_data.size();
      std::cout << " and max_angular_l=";
      std::cout << this->ref_data[0]["max_angular_l"] << std::endl;
    }
    math::SphericalHarmonics harmonics_calculator{};
    size_t max_angular_l = this->ref_data[0]["max_angular_l"];
    harmonics_calculator.precompute(max_angular_l);
    // move precompute out of the function
    for (auto & data : this->ref_data) {
      if (verbose) {
        std::cout << ">> Start loading parameters..." << std::endl;
      }
      // json apparently cannot directly convert to Eigen structures
      std::vector<double> unit_vector_tmp = data["unit_vector"];
      Eigen::Vector3d unit_vector(unit_vector_tmp.data());

      // json apparently cannot directly convert to Eigen structures
      auto harmonics_tmp = data["harmonics"].get<std::vector<double>>();
      math::Vector_t harmonics_ref = Eigen::Map<math::Vector_t>(
          harmonics_tmp.data(), harmonics_tmp.size());

      if (verbose) {
        std::cout << ">> max_angular_l: " << max_angular_l << std::endl;
        std::cout << ">> unit_vector: "
                  << Eigen::Map<Eigen::RowVector3d>(unit_vector.data())
                  << std::endl;
        std::cout << ">> harmonics_ref: " << std::endl;
        std::cout << ">> " << harmonics_ref << std::endl;
      }
      if (verbose) {
        std::cout << ">> Parameters successful loaded." << std::endl;
      }

      // computing with ref parameters with our spherical harmonic
      // harmonics_calculatortion
      harmonics_calculator.calc(unit_vector);
      if (verbose) {
        std::cout << ">> computed harmonics: " << std::endl;
        std::cout << ">> " << harmonics_calculator.get_harmonics() << std::endl;
      }
      double rel_error{
          (harmonics_calculator.get_harmonics() - harmonics_ref).norm()};
      BOOST_CHECK_LE(rel_error, 2 * math::dbl_ftol);
      // TODO(alex) check with main_test_suite -l success and remove
      if (verbose) {
        std::cout << ">> Boost check perfomed." << std::endl;
        std::cout << std::endl;
      }
    }

    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test finished."
                << std::endl;
    }
  }

  /* Rescued from `test_math_utils.cc`: Keep or refer to the script used to
   * generate reference values
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
   */
  BOOST_FIXTURE_TEST_CASE(math_associated_legendre_polynomial_test,
                          SphericalHarmonicsClassRefFixture) {
    verbose = false;
    info = false;
    math::SphericalHarmonics harmonics_calculator{};
    if (info) {
      std::cout << ">> Test math_associated_legendre_polynomial_test "
                << std::endl;
      std::cout << ">> with number_of_unit_vectors=" << this->ref_data.size();
      std::cout << " and max_angular_l=";
      std::cout << this->ref_data[0]["max_angular_l"] << std::endl;
    }
    for (auto & data : this->ref_data) {
      if (verbose) {
        std::cout << ">> Start loading parameters..." << std::endl;
      }
      size_t max_angular_l = data["max_angular_l"];
      // json apparently cannot directly convert to Eigen structures
      std::vector<double> unit_vector_tmp = data["unit_vector"];
      Eigen::Vector3d unit_vector(unit_vector_tmp.data());

      // the matrix is first loaded as vector and then reshaped
      std::vector<double> alps_ref_tmp = data["alps"];
      math::Matrix_t alps_ref = Eigen::Map<math::Matrix_t>(
          alps_ref_tmp.data(), max_angular_l + 1, max_angular_l + 1);

      if (verbose) {
        std::cout << ">> max_angular_l: " << max_angular_l << std::endl;
        std::cout << ">> unit_vector: "
                  << Eigen::Map<Eigen::RowVector3d>(unit_vector.data())
                  << std::endl;
        std::cout << ">> alps_ref: " << std::endl;
        std::cout << alps_ref << std::endl;
      }
      if (verbose) {
        std::cout << ">> Parameters successful loaded." << std::endl;
      }

      // computing with ref parameters with our spherical harmonic
      // harmonics_calculatortion
      harmonics_calculator.precompute(max_angular_l);
      harmonics_calculator.calc(unit_vector);
      if (verbose) {
        std::cout << ">> computed associated legendre polynomial: "
                  << std::endl;
        std::cout << harmonics_calculator.get_assoc_legendre_polynom()
                  << std::endl;
      }
      double rel_error{
          (harmonics_calculator.get_assoc_legendre_polynom() - alps_ref)
              .norm()};
      BOOST_CHECK_LE(rel_error, 10 * math::dbl_ftol);
      // Checks if the additional last column of the associated legendre
      // polynomial matrix contains only zero entries
      BOOST_CHECK_EQUAL(harmonics_calculator.get_assoc_legendre_polynom_raw()
                            .col(max_angular_l + 1)
                            .norm(),
                        0);
    }
    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test finished."
                << std::endl;
    }
  }

  /* TODO(alex) port this test, rescued from test_math_utils.cc, for the new
   * class
   *
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
  */

  BOOST_AUTO_TEST_CASE(spherical_harmonics_gradient_test) {
    // (max) what?! how does this even remain numerically stable?
    // it's total overkill in any case, 10 or even 3 would suffice
    constexpr size_t test_max_angular = 30;
    SphericalHarmonicsGradientsCalculator<test_max_angular>
        harmonics_grad_calc{};
    harmonics_grad_calc.precompute();
    GradientTestFixture fix{
        "reference_data/spherical_harmonics_gradient_test.json"};
    test_gradients(harmonics_grad_calc, fix);
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
