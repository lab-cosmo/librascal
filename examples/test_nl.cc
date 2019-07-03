/**
 * file   test_nl.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2018
 *
 * @brief Example for Neighbour list
 *
 * Copyright  2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/make_structure_manager.hh"
#include "rascal_utility.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap_invariant.hh"
#include "representations/feature_manager_dense.hh"
#include "representations/feature_manager_block_sparse.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

using Representation_t = RepresentationManagerSOAPInvariant<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

int main() {
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  double cutoff{3.};
  json hypers{{"max_radial", 2},
              {"max_angular", 2},
            //   {"soap_type", "PowerSpectrum"},
              {"soap_type", "BiSpectrum"},
              {"inversion_symmetry", false},
              {"normalize", true}};

  json fc_hypers{{"type", "Cosine"},
                 {"cutoff", {{"value", cutoff}, {"unit", "A"}}},
                 {"smooth_width", {{"value", 0.}, {"unit", "A"}}}};
  json sigma_hypers{{"type", "Constant"},
                    {"gaussian_sigma", {{"value", 0.4}, {"unit", "A"}}}};

  hypers["cutoff_function"] = fc_hypers;
  hypers["gaussian_density"] = sigma_hypers;
  hypers["radial_contribution"] = {{"type", "GTO"}};

  json structure{{"filename", filename}};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments",
            {{"cutoff", cutoff}, {"consider_ghost_neighbours", false}}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);
  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>(
          structure, adaptors);
  Representation_t representation{manager, hypers};
  representation.compute();

  size_t inner_size{representation.get_feature_size()};
  FeatureManagerBlockSparse<double> feature{inner_size, hypers};

  feature.push_back(representation);
  auto X{feature.get_feature_matrix_dense()};
  std::cout << "sadfasd" << std::endl;
  auto n_center{feature.sample_size()};
  auto norms = X.colwise().norm();
  std::cout << norms.size() << std::endl;
  for (int icenter{0}; icenter < n_center; icenter++) {
    std::cout << norms[icenter] << std::endl;
  }

  std::cout << "Sample SOAP elements \n"
            << X(0, 0) << " " << X(0, 1) << " " << X(0, 2) << "\n"
            << X(1, 0) << " " << X(1, 1) << " " << X(1, 2) << "\n"
            << X(2, 0) << " " << X(2, 1) << " " << X(2, 2) << "\n";

  // auto& soap_vectors = representation.soap_vectors;

  // for (auto center : manager) {
  //   auto& soap_vector = soap_vectors[center];
  //   for (auto el : soap_vector) {
  //     auto& key = el.first;
  //     auto& soap_vector_by_species = el.second;
  //     std::cout << soap_vector_by_species << std::endl;
  //     break;
  //   }
  // }
  
  auto kernel1 = X.transpose() * X;

  auto kernel2 = dot(feature, feature);

  auto kernel3 = dot(feature);

  auto max1{kernel1.mean()};
  auto max2{kernel2.mean()};
  auto diff{(kernel1 - kernel2).array().abs().matrix().mean()};

  std::cout << max1 << ", " << max2 << ", " << diff << std::endl;

  diff = (kernel2 - kernel3).array().abs().matrix().mean();

  std::cout << diff << std::endl;

  return (0);
}
