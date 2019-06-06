
#include "tests.hh"
#include "test_math.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(MathSphericalHarmonicsTests);
  
  BOOST_FIXTURE_TEST_CASE(math_spherical_harmonics_test, SphericalHarmonicsRefFixture) {
    math::SphericalHarmonics func{};
    for (auto & data : this->ref_data) {
      if (verbose) {
        std::cout << "Start loading parameters." << std::endl;
      }
      size_t max_angular_l = data["max_angular_l"] ;
      // the reader is invited to do these two lines in one line, if the reader
      // is able to
      std::vector<double> unit_vector = data["unit_vector"];
      Eigen::Vector3d curr_unit_vector(unit_vector.data());
      
      std::vector<double> curr_harmonics = data["harmonics"];
      math::Vector_t curr_harmonics_ref =
        Eigen::Map<math::Vector_t>(
            curr_harmonics.data(), curr_harmonics.size());
      if (verbose) {
        std::cout << "max_angular_l: " << max_angular_l << std::endl;
        std::cout << "unit_vector: " << //curr_unit_vector << std::endl;
            Eigen::Map<Eigen::RowVector3d>(curr_unit_vector.data()) << std::endl;
        std::cout << "harmonics_ref: " << curr_harmonics_ref << std::endl;
        //    Eigen::Map<Eigen::Vector3d>(curr_unit_vector.data()) << std::endl;
        //std::cout << "harmonics: " << 
        //    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        //    curr_harmonics.data(), curr_harmonics.size()) << std::endl;
        std::cout << "Parameters successful loaded" << std::endl;
      }

      // computing with ref parameters with our spherical harmonic function
      
      func.precompute(max_angular_l);
      func.compute(curr_unit_vector);
      if (verbose) {        
        std::cout << "harmonics: " << func.get_harmonics() << std::endl;
      }
      //double rel_error{(curr_harmonics_ref -  func.get_harmonics()).norm()};
      auto && harmonics_tmp{math::compute_spherical_harmonics(curr_unit_vector, max_angular_l)};
      //if (verbose) {        
      //  std::cout << "harmonics_tmp: " << harmonics_tmp << std::endl;
      //    std::cout << harmonics_tmp.isRowMajor << std::endl; 
      //  //for (size_t i{0}; i<curr_harmonics_ref.size(); i++) {          
      //  //  std::cout << curr_harmonics_ref(i) << " "<<  harmonics_tmp(i) << std::endl;
      //  //  //double tmp = (curr_harmonics_ref(i) - harmonics_tmp(i));
      //  //  //std::cout << "error" << tmp << std::endl;
      //  //}
      //}
      double rel_error{(func.get_harmonics() -  harmonics_tmp).norm()};
      //double rel_error{(func.get_harmonics() -  curr_harmonics_ref).norm()};
      BOOST_CHECK_LE(rel_error, 15 * math::dbl_ftol);
      //double rel_error{(curr_harmonics_ref -  harmonics).norm()};
      //BOOST_CHECK_LE(rel_error, 15 * math::dbl_ftol);
    }
    // does not work
    //const auto & max_angular_l{this->ref_data.at("max_angular_l").template get<double>()};
    
    //for (size_t i{0}; i< max_angular_l.capacity(); i++) {
    //  // loading ref parameters and results
    //  auto && curr_max_angular_l = max_angular_l.at(i);
    //  if (verbose) {
    //    std::cout << "Start loading parameters" << std::endl;
    //  }
    //  Eigen::Vector3d curr_unit_vector(unit_vector.at(i).data());
    //  //size_t harmonics_size{(curr_max_angular_l+1)*(curr_max_angular_l+1)};
    //  auto && curr_harmonics = harmonics.at(i);
    //  math::Vector_t curr_harmonics_ref =
    //    Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
    //        curr_harmonics.data(), curr_harmonics.size());

    //  // computing with ref parameters with our spherical harmonic function
    //  func.precompute(curr_max_angular_l);
    //  func.compute(curr_unit_vector);
    //  double rel_error{(curr_harmonics_ref -  func.get_harmonics()).norm()};

    //  BOOST_CHECK_LE(rel_error, 15 * math::dbl_ftol);
    //}

  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
