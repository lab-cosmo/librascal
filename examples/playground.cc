/**
 * file   playground.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2018
 *
 * @brief File to build and test new functionality
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
#include "structure_managers/adaptor_center_contribution.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/make_structure_manager.hh"
#include "structure_managers/structure_manager_collection.hh"

#include "rascal_utility.hh"
#include "representations/calculator_sorted_coulomb.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"

#include "models/kernels.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>
#include <random>

using namespace rascal;  // NOLINT

// template <typename Manager, template <class> class...
// AdaptorImplementationPack> struct Test {
//   using ManagerTypeHolder_t =
//       StructureManagerTypeHolder<StructureManagerCenters,
//       AdaptorNeighbourList,
//                                  AdaptorStrict>;
//   using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
//   void operator()() {
//     std::cout << internal::GetTypeName<ManagerTypeList_t>() << std::endl;
//   }
// };

// // using Representation_t = CalculatorSphericalInvariants;
// using ManagerTypeHolder_t =
//     StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
//                                AdaptorCenterContribution, AdaptorStrict>;
// using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
// using Manager_t = typename ManagerTypeHolder_t::type;
// using ManagerCollection_t =
//     typename TypeHolderInjector<ManagerCollection, ManagerTypeList_t>::type;
// // using Manager_t =
// // AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>;
// using Representation_t = CalculatorSphericalInvariants;
// using Property_t = typename Representation_t::template Property_t<Manager_t>;
// // using ManagerCollection_t = ManagerCollection<>;
// using Test1 = typename TypeHolderInjector<Test, ManagerTypeList_t>::type;

// template <typename T, size_t Order, int NbRow = 1, int NbCol = 1>
// using Prop_t = Property<T, Order, 1, Manager_t, NbRow, NbCol>;

int main() {
  // Test1()();
  // std::string filename{"reference_data/dft-smiles_500.ubjson"};
  // std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  std::string filename{"reference_data/methane.json"};
  // std::string filename{"reference_data/diamond_cubic.json"};
  std::string rep_id{"pp"};

  double cutoff{3.};
  // json hypers{{"max_radial", 6},
  //             {"max_angular", 6},
  //             {"soap_type", "PowerSpectrum"},
  //             {"normalize", true},
  //             {"identifier",rep_id}};
  json hypers{{"max_radial", 8},
              {"max_angular", 6},
              {"soap_type", "PowerSpectrum"},
              // {"soap_type", "BiSpectrum"},
              {"inversion_symmetry", true},
              {"normalize", true},
              {"compute_gradients", false}};

  json fc_hypers{{"type", "Cosine"},
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
           {"initialization_arguments",
            {{"cutoff", cutoff},
             {"consider_ghost_neighbours", false},
             {"skin", 0.}}}};
  json ad1b{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad1b);
  adaptors.emplace_back(ad2);
  // ManagerCollection_t collection{adaptors};
  // collection.add_structures(filename, 0, 10);
  // std::cout << collection.size() << std::endl;

  // Representation_t representation{hypers};
  // representation.compute(collection);

  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structure, adaptors);
  // auto manager =
  //     make_structure_manager_stack<StructureManagerCenters,
  //                                  AdaptorNeighbourList,
  //                                  AdaptorStrict>(
  //         structure, adaptors);
  std::cout << "n_centers: " << manager->size() << std::endl;
  for (auto center : manager) {
    auto ctag = center.get_atom_tag();
    std::cout << "Center: " << ctag << std::endl;
    auto atom_ii = center.get_atom_ii();
    auto atom_ii_tag = atom_ii.get_atom_tag_list();
    auto atom_ii_ids = atom_ii.get_cluster_indices();

    for (auto neigh : center) {
      auto tag_list = neigh.get_atom_tag_list();
      auto dist = manager->get_distance(neigh);

      auto atom_j = neigh.get_atom_j();
      auto atom_j_tag = atom_j.get_atom_tag_list();
      auto atom_j_ids = atom_j.get_cluster_indices();
      std::cout << "neigh: " << tag_list[0] << ", " << tag_list[1] << ", "
                << dist << " tag_ii: " << atom_ii_tag[0] << ", "
                << atom_ii_tag[1] << ", " << atom_ii_ids[0]
                << " tag_j: " << atom_j_tag[0] << ", " << atom_j_ids[0]
                << std::endl;
    }
  }

  // // using Prop_t = typename
  // //
  // StructureManager<AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>::template
  // // Property_t<int, 1, 1>;

  // auto prop = Prop_t<int, 1>(*manager);
  // int ii{0};
  // prop.resize();
  // for (auto center : manager) {
  //   prop[center] = ii;
  //   ++ii;
  // }

  // for (auto center : manager) {
  //   auto && center_tag = center.get_atom_tag();
  //   std::cout << "Center prop: " << center_tag << " -- " << prop[center]
  //             << " -- " << manager->get_position(center_tag).transpose()
  //             << std::endl;
  //   for (auto neigh : center) {
  //     std::cout << "typeid(neigh) " << typeid(neigh).name() << std::endl;
  //     auto atom_j = neigh.get_atom_j();
  //     auto && atom_j_tag = atom_j.get_atom_tag();
  //     std::cout << "atom_j prop: " << atom_j_tag << " -- " << prop[atom_j]
  //               << " -- " << manager->get_position(atom_j_tag).transpose()
  //               << std::endl;
  //   }
  // }

  // for (auto center : manager) {
  //   auto&& center_tag = center.get_atom_tag();
  //   std::cout << "Center tag: "<< center_tag << " -- " <<
  //   manager->get_atom_index(center_tag) << " -- " << center.get_atom_type()
  //   << " -- " << center.get_position().transpose() <<  std::endl; for (auto
  //   neigh : center) {
  //     auto&& neigh_tag = neigh.get_atom_tag();
  //     std::cout << "Neigh tag: "<< neigh_tag << " -- " <<
  //     manager->get_atom_index(neigh_tag) << " -- " << neigh.get_atom_type()
  //     << " -- " << neigh.get_position().transpose() <<  std::endl; auto
  //     atom_j = neigh.get_atom_j(); auto&& atom_j_tag = atom_j.get_atom_tag();
  //     std::cout << "atom_j tag: "<< atom_j_tag << " -- " <<
  //     manager->get_atom_index(atom_j_tag) << " -- " <<
  //     manager->get_atom_type(atom_j_tag) << " -- " <<
  //     manager->get_position(atom_j_tag).transpose() <<  std::endl;
  //   }
  // }

  // ManagerCollection_t collectionA{adaptors};
  // collectionA.add_structures(filename, 0, 10);
  // std::cout << collectionA.size() << std::endl;
  // ManagerCollection_t collectionB{adaptors};
  // collectionB.add_structures(filename, 0, 10);
  // std::cout << collectionB.size() << std::endl;

  // Representation_t representation{hypers};
  // representation.compute(collectionA);
  // representation.compute(collectionB);

  // json kernel_hypers{{"zeta", 2},
  //                     {"target_type", "Structure"},
  //                     {"name", "Cosine"}};
  // Kernel kernel{kernel_hypers};
  // // for (auto& manager : collectionA) {
  // //   for (auto center : manager) {
  // //     std::cout << center.get_position().transpose() << std::endl;
  // //   }
  // //   std::cout  << std::endl;
  // // }
  // auto mat = kernel.compute(representation, collectionA, collectionB);
  // std::cout << mat << std::endl;

  // // json kernel_hypers_local{{"zeta", 2},
  // //                     {"target_type", "Atom"}};
  // // Kernel<internal::KernelType::Cosine> kernel_local{kernel_hypers_local};

  // // auto mat_local = kernel_local.compute(representation, collectionA,
  // collectionA);
  // // std::cout << mat_local << std::endl;

  // // auto property_name{representation.get_name()};
  // // auto&& property{manager->template
  // get_validated_property_ref<Property_t>(property_name)};

  // // auto test_representation{property.get_dense_rep()};

  // TODELETE
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

  // auto test_representation{property.get_dense_feature_matrix()};

  return (0);
}
