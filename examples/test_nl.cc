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

struct ordering {
  template <typename T>
  bool operator()(std::pair<size_t, T> const & a,
                  std::pair<size_t, T> const & b) {
    return *(a.second) < *(b.second);
  }
};

template <typename T>
auto sort_with_ordering(T & Index) {
  using myiter = typename T::const_iterator;
  using ret_t = std::vector<std::pair<size_t, myiter>>;
  ret_t order(Index.size());

  size_t n = 0;
  for (myiter it = Index.begin(); it != Index.end(); ++it, ++n) {
    order[n] = std::make_pair(n, it);
  }

  std::sort(order.begin(), order.end(), ordering());

  return order;
}

template <typename T, typename V>
std::vector<T> sort_from_ref(const std::vector<T> & in, const V & order) {
  std::vector<T> ret(in.size());

  for (size_t i = 0; i < in.size(); ++i) {
    ret[i] = in[order[i].first];
  }

  return ret;
}

int main() {
  std::string filename{"reference_data/small_molecule.json"};
  double cutoff{3.};
  json hypers{{"max_radial", 6},
              {"max_angular", 6},
              {"soap_type", "PowerSpectrum"},
              {"normalize", true},
              {"compute_gradients", true}};

  json fc_hypers{{"type", "Cosine"},
                 {"cutoff", {{"value", cutoff}, {"unit", "AA"}}},
                 {"smooth_width", {{"value", 0.}, {"unit", "AA"}}}};
  json sigma_hypers{{"type", "Constant"},
                    {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};

  hypers["cutoff_function"] = fc_hypers;
  hypers["gaussian_density"] = sigma_hypers;
  hypers["radial_contribution"] = {{"type", "GTO"}};

  // json hypers{{"central_decay", 0.5},
  //             {"interaction_cutoff", 10.},
  //             {"interaction_decay", 0.5},
  //             {"size", 14},
  //             {"sorting_algorithm", "distance"}};

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

  Representation_t representation{manager, hypers};
  representation.compute();

  auto rep_full = representation.get_representation_full();

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
  std::vector<std::vector<double>> types_ref{};
  std::vector<std::vector<double>> types{};

  size_t i_center{0};
  auto mm0 = extract_underlying_manager<0>(manager);
  for (auto center : manager) {
    if (is_center_atom(i_center)) {
      distances_ref.emplace_back();
      types_ref.emplace_back();
      std::cout << "center_atom: " << center.get_atom_tag() << " -- "
                << mm0->get_atom_index(center) << " -- "
                << center.get_position().transpose() << std::endl;
      for (auto neigh : center) {
        auto neigh_tag = neigh.get_atom_tag();
        auto neigh_type = neigh.get_atom_type();
        std::cout << "neigh_atom: " << neigh_type << " -- "
                  << manager->get_neighbour_atom_tag(center, neigh.get_index())
                  << " -- " << neigh.get_position().transpose() << " -- "
                  << manager->get_position(neigh_tag).transpose() << std::endl;
        auto dist{(neigh.get_position() - center.get_position()).norm()};
        distances_ref.back().push_back(dist);
        types_ref.back().push_back(neigh_type);
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
    types.emplace_back();
    std::cout << "center_atom: " << center.get_atom_tag() << " -- "
              << mm1->get_atom_index(center) << " -- "
              << center.get_position().transpose() << std::endl;
    for (auto neigh : center) {
      auto neigh_tag = neigh.get_atom_tag();
      auto neigh_type = neigh.get_atom_type();
      neigh.get_position().transpose();
      std::cout << "neigh_atom: " << neigh_type << " -- "
                << manager->get_neighbour_atom_tag(center, neigh.get_index())
                << " -- " << neigh.get_position().transpose() << " -- "
                << manager->get_position(neigh_tag).transpose() << std::endl;
      auto dist{(neigh.get_position() - center.get_position()).norm()};
      distances.back().push_back(dist);
      types.back().push_back(neigh_type);
      // distances.push_back(manager_no_center->get_distance(neigh));
    }
  }

  std::cout << "is_center_atom: " << is_center_atom.transpose() << std::endl;
  i_center = 0;
  for (; i_center < manager->size(); ++i_center) {
    auto ref_order = sort_with_ordering(distances_ref[i_center]);
    auto order = sort_with_ordering(distances[i_center]);

    auto distance_ref = sort_from_ref(distances_ref[i_center], ref_order);
    auto distance = sort_from_ref(distances[i_center], order);

    auto type_ref = sort_from_ref(types_ref[i_center], ref_order);
    auto type = sort_from_ref(types[i_center], order);

    std::cout << "Center: " << i_center << std::endl;
    std::cout << "sizes: " << distances_ref[i_center].size() << ", "
              << distances[i_center].size() << std::endl;
    for (size_t i_d{0}; i_d < distances[i_center].size(); i_d++) {
      std::cout << type_ref[i_d] << "\t" << type[i_d] << "\t"
                << std::abs(distance_ref[i_d] - distance[i_d]) << "\t"
                << distance_ref[i_d] << "\t" << distance[i_d] << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << " ###################REP################## " << std::endl;

  Representation_t representation_no_center{manager, hypers};
  representation_no_center.compute();

  auto rep_no_center = representation_no_center.get_representation_full();

  size_t i_no_center{0};
  std::cout << manager->size();
  for (i_center = 0; i_center < n_atoms; ++i_center) {
    std::cout << "Center idx: " << i_center << std::endl;
    if (is_center_atom(i_center)) {
      auto row_full = rep_full.col(i_center).eval();
      auto row_no_center = rep_no_center.col(i_no_center).eval();
      auto diff = (row_full - row_no_center).norm();
      std::cout << "Center idx: " << i_center << " Diff: " << diff << std::endl;
      i_no_center++;
    }
  }


  size_t inner_size{representation.get_feature_size()};
  FeatureManagerBlockSparse<double> feature{inner_size, hypers};

  feature.push_back(representation);
  auto X{feature.get_feature_matrix_dense()};
  std::cout << "sadfasd" << std::endl;
  // auto n_center{feature.sample_size()};
  auto norms = X.colwise().norm();
  std::cout << norms.size() << std::endl;
  // for (int icenter{0}; icenter < n_center; icenter++) {
  //   std::cout << norms[icenter] << std::endl;
  // }

  auto kernel1 = X.transpose() * X;

  auto kernel2 = dot(feature, feature);

  auto kernel3 = dot(feature);

  // auto max1{kernel1.mean()};
  // auto max2{kernel2.mean()};
  // auto diff{(kernel1 - kernel2).array().abs().matrix().mean()};
  std::cout << kernel1.mean() << ", " << kernel1.minCoeff() << ", "
            << kernel1.maxCoeff() << std::endl;
  // std::cout << max1 << ", " << max2 << ", " << diff << std::endl;

  // diff = (kernel2 - kernel3).array().abs().matrix().mean();

  // std::cout << diff << std::endl;

  return (0);
}
