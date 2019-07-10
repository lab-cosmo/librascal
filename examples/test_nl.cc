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
#include "representations/representation_manager_soap.hh"
#include "representations/feature_manager_dense.hh"
#include "representations/feature_manager_block_sparse.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>
#include <random>

// using namespace std;
using namespace rascal;  // NOLINT

using Representation_t = RepresentationManagerSOAP<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;
using ArrayB_t = AtomicStructure<3>::ArrayB_t;

int main() {
  std::string filename{"reference_data/small_molecule.json"};
  double cutoff{3.};
  json hypers{{"max_radial", 6},
              {"max_angular", 6},
              {"soap_type", "PowerSpectrum"},
              {"normalize", true}};

  json fc_hypers{{"type", "Cosine"},
                 {"cutoff", {{"value", cutoff}, {"unit", "AA"}}},
                 {"smooth_width", {{"value", 0.}, {"unit", "AA"}}}};
  json sigma_hypers{{"type", "Constant"},
                    {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};

  hypers["cutoff_function"] = fc_hypers;
  hypers["gaussian_density"] = sigma_hypers;
  hypers["radial_contribution"] = {{"type", "GTO"}};

  json structure{};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments",
            {{"cutoff", cutoff}, {"consider_ghost_neighbours", true}}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);

  AtomicStructure<3> atomic_structure{};
  atomic_structure.set_structure(filename);
  structure = atomic_structure;

  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>(
          structure, adaptors);

  // manager->update(atomic_structure);
  std::mt19937_64 rng{1242484542};

  auto n_atoms = atomic_structure.get_number_of_atoms();
  std::uniform_int_distribution<int> uni(1, n_atoms - 2);
  // auto n_flips = uni(rng);
  // int n_flips{4};
  // for (int i_it{0}; i_it < n_flips; ++i_it) {
  //   auto i_idx = uni(rng);
  //   atomic_structure.is_a_center_atom(i_idx) = false;
  // }
  // atomic_structure.is_a_center_atom = false;
  // // atomic_structure.is_a_center_atom(0) = true;
  // atomic_structure.is_a_center_atom(4) = true;
  atomic_structure.is_a_center_atom(2) = false;
  atomic_structure.is_a_center_atom(3) = false;

  atomic_structure.is_a_center_atom(6) = false;
  atomic_structure.is_a_center_atom(7) = false;
  // atomic_structure.is_a_center_atom(10) = false;
  // atomic_structure.is_a_center_atom(12) = false;
  // atomic_structure.is_a_center_atom(13) = false;
  // atomic_structure.is_a_center_atom(14) = false;
  // atomic_structure.is_a_center_atom(15) = false;
  // atomic_structure.is_a_center_atom(16) = false;

  auto & is_center_atom = atomic_structure.is_a_center_atom;

  std::vector<std::vector<double>> distances_ref{};
  std::vector<std::vector<double>> distances{};

  size_t i_center{0};
  auto mm0 = extract_underlying_manager<0>(manager);
  for (auto center : manager) {
    if (is_center_atom(i_center)) {
      distances_ref.emplace_back();
      std::cout << "center_atom: " << center.get_atom_tag() << " -- "
                << mm0->get_atom_index(center) << " -- "
                << center.get_position().transpose() << std::endl;
      for (auto neigh : center) {
        auto neigh_tag = neigh.get_atom_tag();
        std::cout << "neigh_atom: " << neigh_tag << " -- "
                  << manager->get_neighbour_atom_tag(center, neigh.get_index())
                  << " -- " << neigh.get_position().transpose() << " -- "
                  << manager->get_position(neigh_tag).transpose() << std::endl;
        auto dist{(neigh.get_position() - center.get_position()).norm()};
        distances_ref.back().push_back(dist);
        // distances_ref.push_back(manager->get_distance(neigh));
      }
    }
    i_center++;
  }
  std::cout << std::endl;
  manager->update(atomic_structure);
  // std::vector<int> ids{{191,124,44,127}};
  // for (auto& idx : ids) {
  //   std::cout << "idx: "<< idx << " -- "<<
  //   manager->get_position(idx).transpose() << std::endl;
  // }
  // std::cout << std::endl;

  // for (int idx{0}; idx < ghost_pos.cols(); ++idx) {
  //   std::cout << "idx: "<< idx << " -- "<< ghost_pos.col(idx).transpose() <<
  //   std::endl;
  // }
  // std::cout << std::endl;

  auto mm1 = extract_underlying_manager<0>(manager);
  for (auto center : manager) {
    distances.emplace_back();
    std::cout << "center_atom: " << center.get_atom_tag() << " -- "
              << mm1->get_atom_index(center) << " -- "
              << center.get_position().transpose() << std::endl;
    for (auto neigh : center) {
      auto neigh_tag = neigh.get_atom_tag();
      neigh.get_position().transpose();
      std::cout << "neigh_atom: " << neigh_tag << " -- "
                << manager->get_neighbour_atom_tag(center, neigh.get_index())
                << " -- " << neigh.get_position().transpose() << " -- "
                << manager->get_position(neigh_tag).transpose() << std::endl;
      auto dist{(neigh.get_position() - center.get_position()).norm()};
      distances.back().push_back(dist);
      // distances.push_back(manager_no_center->get_distance(neigh));
    }
  }

  std::cout << "is_center_atom: " << is_center_atom.transpose() << std::endl;
  i_center = 0;
  for (; i_center < manager->size(); ++i_center) {
    std::sort(distances_ref[i_center].begin(), distances_ref[i_center].end());
    std::sort(distances[i_center].begin(), distances[i_center].end());
    std::cout << "Center: " << i_center << std::endl;
    std::cout << "sizes: " << distances_ref[i_center].size() << ", "
              << distances[i_center].size() << std::endl;
    for (size_t i_d{0}; i_d < distances[i_center].size(); i_d++) {
      std::cout << std::abs(distances_ref[i_center][i_d] -
                            distances[i_center][i_d])
                << "\t" << distances_ref[i_center][i_d] << "\t"
                << distances[i_center][i_d] << std::endl;
    }
  }

  // AtomicStructure<3> atomic_structure2{atomic_structure};

  // atomic_structure2.positions(0, 0) += 0.5;

  // manager->update(atomic_structure2);

  // Representation_t representation{manager, hypers};
  // representation.compute();

  // size_t inner_size{representation.get_feature_size()};
  // FeatureManagerBlockSparse<double> feature{inner_size, hypers};

  // feature.push_back(representation);
  // auto X{feature.get_feature_matrix_dense()};
  // std::cout << "sadfasd" << std::endl;
  // auto n_center{feature.sample_size()};
  // auto norms = X.colwise().norm();
  // std::cout << norms.size() << std::endl;
  // for (int icenter{0}; icenter < n_center; icenter++) {
  //   std::cout << norms[icenter] << std::endl;
  // }

  // auto kernel1 = X.transpose() * X;

  // auto kernel2 = dot(feature, feature);

  // auto kernel3 = dot(feature);

  // auto max1{kernel1.mean()};
  // auto max2{kernel2.mean()};
  // auto diff{(kernel1 - kernel2).array().abs().matrix().mean()};

  // std::cout << max1 << ", " << max2 << ", " << diff << std::endl;

  // diff = (kernel2 - kernel3).array().abs().matrix().mean();

  // std::cout << diff << std::endl;

  return (0);
}
