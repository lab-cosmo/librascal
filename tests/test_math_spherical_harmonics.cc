
#include "tests.hh"
#include "test_math.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(MathSphericalHarmonicsTests);

  // We assume that max_angular_l does not change in the reference data.
  // TODO(alex) make this more clear by changing the data structure of the
  // reference data
  // TODO(alex) replace info and verbose usage of VerbosityValue NORMAL DEBUG
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
      std::vector<double> harmonics_tmp = data["harmonics"];
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

      // computing with ref parameters with our spherical harmonic harmonics_calculatortion
      harmonics_calculator.calc(unit_vector);
      if (verbose) {
        std::cout << ">> computed harmonics: " << std::endl;
        std::cout << ">> " << harmonics_calculator.get_harmonics() << std::endl;
      }
      double rel_error{(harmonics_calculator.get_harmonics() - harmonics_ref).norm()};
      BOOST_CHECK_LE(rel_error, 2*math::dbl_ftol);
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

  BOOST_FIXTURE_TEST_CASE(math_associated_legendre_polynomial_test,
                          SphericalHarmonicsClassRefFixture) {
    verbose = false;
    info = false;
    math::SphericalHarmonics harmonics_calculator{};
    if (info) {
      std::cout << ">> Test math_associated_legendre_polynomial_test " << std::endl;
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

      // computing with ref parameters with our spherical harmonic harmonics_calculatortion
      harmonics_calculator.precompute(max_angular_l);
      harmonics_calculator.calc(unit_vector);
      if (verbose) {
        std::cout << ">> computed associated legendre polynomial: "
                  << std::endl;
        std::cout << harmonics_calculator.get_assoc_legendre_polynom() << std::endl;
      }
      double rel_error{(harmonics_calculator.get_assoc_legendre_polynom() - alps_ref).norm()};
      BOOST_CHECK_LE(rel_error, 10 * math::dbl_ftol);
      // Checks if the additional last column of the associated legendre polynomial matrix contains only zero entries
      BOOST_CHECK_EQUAL(harmonics_calculator.get_assoc_legendre_polynom_raw().col(max_angular_l+1).norm(), 0);
      
    }
    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test finished."
                << std::endl;
    }
  }

  BOOST_AUTO_TEST_CASE(spherical_harmonics_gradient_test) {
    constexpr size_t test_max_angular = 30;
    SphericalHarmonicsClassFunctionHolder<test_max_angular> harmonics_calculator{};
    harmonics_calculator.precompute();
    GradientTestFixture fix{
        "reference_data/spherical_harmonics_gradient_test.json"};
    test_gradients(harmonics_calculator, fix);
  }


  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
