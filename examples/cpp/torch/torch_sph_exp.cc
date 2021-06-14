/**
 * @file   examples/spherical_expansion_example.cc
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   26 June 2019
 *
 * @brief  Example for computing the spherical expansion
 *
 * Copyright © 2018 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/utils/utils.hh"
#include "rascal/math/utils.hh"

#include <torch/torch.h>
#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <string>
#include <stdexcept>

using namespace rascal;  // NOLINT

using Matrix_t = math::Matrix_t;

using Representation_t = CalculatorSphericalExpansion;
using Manager_t = AdaptorStrict<
    AdaptorCenterContribution<AdaptorNeighbourList<StructureManagerCenters>>>;
using Prop_t = typename CalculatorSphericalExpansion::Property_t<Manager_t>;
using PropGrad_t =
    typename CalculatorSphericalExpansion::PropertyGradient_t<Manager_t>;


// adapted from https://discuss.pytorch.org/t/data-transfer-between-libtorch-c-and-eigen/54156/6
torch::Tensor eigenMatrixToTorchTensor(Matrix_t eigen_mat) {
    // one can play more around with the options, not sure at the moment
    auto options =
      torch::TensorOptions()
        .dtype(torch::kDouble)
        .layout(torch::kStrided)
        .device(torch::kCPU) // kCUDA, 1
        .requires_grad(true);
     
    auto tensor = torch::empty({eigen_mat.rows(), eigen_mat.cols()}, options);
    double* data = tensor.data_ptr<double>();

    // TODO check from_blop 
    Eigen::Map<Matrix_t> eigen_mat_cast(data, tensor.size(0), tensor.size(1));
    eigen_mat_cast = eigen_mat.cast<double>();
    return tensor;
}


int main(int argc, char * argv[]) {
  if (argc < 4) {
    std::cerr << "Must provide {atomic structure json filename} {max_radial} {max_angular}";
    std::cerr << std::endl;
    return -1;
  }

  std::string filename{argv[1]};
  int max_radial{std::atoi(argv[2])};
  int max_angular{std::atoi(argv[3])};

  double cutoff{2.};
  json hypers{
      {"max_radial", max_radial}, {"max_angular", max_angular}, {"compute_gradients", true}};

  json fc_hypers{{"type", "ShiftedCosine"},
                 {"cutoff", {{"value", cutoff}, {"unit", "AA"}}},
                 {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}};
  json sigma_hypers{{"type", "Constant"},
                    {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};

  hypers["cutoff_function"] = fc_hypers;
  hypers["gaussian_density"] = sigma_hypers;
  hypers["radial_contribution"] = {{"type", "GTO"}};

  json structure{{"filename", filename}};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  json ad1b{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad1b);
  adaptors.emplace_back(ad2);

  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structure, adaptors);

  Representation_t representation{hypers};
  representation.compute(manager);

  int nb_centers{0};
  std::set<int> species{};
  for (auto center : manager) {
    auto ii_pair = center.get_atom_ii();
    species.insert(center.get_atom_type());
    nb_centers++;
  }
  std::cout << "nb_centers: " << nb_centers << std::endl;
  std::cout << "nb_species: " << species.size() << std::endl;
  std::cout << "max_radial: " << max_radial << std::endl;
  std::cout << "max_angular: " << max_angular << std::endl;
  auto && expansions_coefficients{
      *manager->template get_property<Prop_t>(representation.get_name())};
  Matrix_t feature_matrix = expansions_coefficients.get_features();
  //std::cout << "nb_features: nb_species*max_radial*(max_angular+1)**2" << std::endl;
  //std::cout << "nb_features: " << feature_matrix.cols() << std::endl;

  torch::Tensor feature_matrix_torch = eigenMatrixToTorchTensor(feature_matrix);

  //std::cout << std::endl;
  std::cout << "EIGEN Matrix (" << feature_matrix.rows() << ", " << feature_matrix.cols() << ")" << std::endl;
  std::cout << feature_matrix << std::endl;
  std::cout << std::endl;
  std::cout << "TORCH Matrix" << std::endl;
  std::cout << feature_matrix_torch  << std::endl;
  
  auto && expansions_coefficients_gradient{
      *manager->template get_property<PropGrad_t>(
          representation.get_gradient_name())};

  Matrix_t feature_gradient_matrix = expansions_coefficients_gradient.get_features_gradient();
  std::cout << std::endl;
  std::cout << "EIGEN Gradient Matrix (" << feature_gradient_matrix.rows() << ", " << feature_gradient_matrix.cols() << ")" << std::endl;
  //std::cout << feature_gradient_matrix << std::endl;
}
