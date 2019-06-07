
#include "tests.hh"
#include "test_math.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(MathSphericalHarmonicsTests);
  
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test, SphericalHarmonicsClassRefFixture) {
    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test " << std::endl;
      std::cout << ">> with number_of_unit_vectors=" << this->ref_data.size();
      std::cout << " and max_angular_l=";
      std::cout << this->ref_data[0]["max_angular_l"] << std::endl;
    }
    math::SphericalHarmonics func{};
    for (auto & data : this->ref_data) {
      if (verbose) {
        std::cout << ">> Start loading parameters..." << std::endl;
      }
      size_t max_angular_l = data["max_angular_l"] ;
      // the reader is invited to do these two lines in one line, if the reader
      // is able to
      std::vector<double> unit_vector_tmp = data["unit_vector"];
      Eigen::Vector3d unit_vector(unit_vector_tmp.data());
      
      // the reader is invited to do these two lines in one line, if the reader
      // is able to
      std::vector<double> harmonics_tmp = data["harmonics"];
      math::Vector_t harmonics_ref =
        Eigen::Map<math::Vector_t>(
            harmonics_tmp.data(), harmonics_tmp.size());

      if (verbose) {
        std::cout << ">> max_angular_l: " << max_angular_l << std::endl;
        std::cout << ">> unit_vector: " << 
            Eigen::Map<Eigen::RowVector3d>(unit_vector.data()) << std::endl;
        std::cout << ">> harmonics_ref: " << std::endl;
        std::cout << ">> " << harmonics_ref << std::endl;
      }
      if (verbose) {
        std::cout << ">> Parameters successful loaded." << std::endl;
      }

      // computing with ref parameters with our spherical harmonic function
      func.precompute(max_angular_l);
      func.compute(unit_vector);
      if (verbose) {        
        std::cout << ">> computed harmonics: " << std::endl;
        std::cout << ">> " << func.get_harmonics() << std::endl;
      }
      double rel_error{(func.get_harmonics() -  harmonics_ref).norm()};
      BOOST_CHECK_LE(rel_error, math::dbl_ftol);
      if (verbose) {
        std::cout << ">> Boost check perfomed." << std::endl;
        std::cout << std::endl;
      }

    }
    if (info) {
      std::cout << ">> Test math_spherical_harmonics_test finished." << std::endl;
    }

  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
