#include <iostream>
#include <fstream>

#include "representations/representation_manager_spherical_expansion.hh"

using namespace rascal;  // NOLINT
using Vector_Ref = internal::RadialContributionBase::Vector_Ref;

int main () {
  int nb_points{100};
  //Vector_t points(11);
  //points << 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.;
  math::Vector_t points{math::Vector_t::LinSpaced(nb_points, 0, 10)};
  double sigma{0.5};

  int max_radial{21};
  int max_angular{max_radial-1};
  json fc_hypers{
       {"type", "Constant"},
       {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
      };
  json hypers{{"gaussian_density", fc_hypers},
            {"max_radial", max_radial},
            {"max_angular", max_angular},
            {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
  };


  math::Matrix_t results(points.size(), max_radial);
  
  auto radial_contr{internal::RadialContribution<internal::RadialBasisType::GTO>(hypers)};
  std::ofstream myfile;
  std::string filename = "radial_contr_logs/radial_contr_nmax"+ std::to_string(max_radial) + "lmax" + std::to_string(max_angular) +".log";
  std::cout << "Saving results into " + filename << std::endl;
  myfile.open(filename);
  for (int i{0}; i<points.size(); i++) {
    // TODO should be aware that we only save part of the matrix here
    myfile << radial_contr.compute_contribution<internal::AtomicSmearingType::Constant>(points(i), sigma) << std::endl;
  }
  myfile.close();
  std::cout << "Saved results into " + filename << std::endl;
  return 0;
}
