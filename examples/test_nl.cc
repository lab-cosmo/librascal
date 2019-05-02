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
#include "representations/feature_manager_dense.hh"
#include "basic_types.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

using Representation_t = RepresentationManagerSphericalExpansion<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

template <typename StructureManagerTypeHolder>
decltype(auto) wrap_factory() {
  return [](const json & a, const json & b) {
    return make_structure_manager_stack_with_hypers_and_typeholder<
        typename StructureManagerTypeHolder::type_list>::apply(a, b);
  };
}

struct TestData {
  using ManagerTypeHolder_t =
      StructureManagerTypeHolder<StructureManagerCenters,
                                  AdaptorNeighbourList, AdaptorStrict>;

  TestData() = default;

  void get_ref(const std::string & ref_filename) {
    std::vector<std::uint8_t> ref_data_ubjson;
    internal::read_binary_file(ref_filename, ref_data_ubjson);
    this->ref_data = json::from_ubjson(ref_data_ubjson);
    auto filenames = this->ref_data.at("filenames").get<std::vector<std::string>>();
    auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

    for (auto && filename : filenames) {
      for (auto && cutoff : cutoffs) {
        json parameters;
        json structure{{"filename", filename}};
        json adaptors;
        json ad1{
            {"name", "AdaptorNeighbourList"},
            {"initialization_arguments",
              {{"cutoff", cutoff},
              {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
        json ad2{{"name", "AdaptorStrict"},
                  {"initialization_arguments", {{"cutoff", cutoff}}}};
        adaptors.emplace_back(ad1);
        adaptors.emplace_back(ad2);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);


      }
    }
  }

  ~TestData() = default;

  const bool consider_ghost_neighbours{false};
  json ref_data{};
  json factory_args{};

};

int main() {

  using ManagerTypeHolder_t =
  StructureManagerTypeHolder<StructureManagerCenters,
          AdaptorNeighbourList, AdaptorStrict>;
  using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
  using Manager_t = typename ManagerTypeHolder_t::type;
  using Std2DArray_t = std::vector<std::vector<double>>;

  auto dd{TestData()};
  std::string filename{"reference_data/sorted_coulomb_reference.ubjson"};
  dd.get_ref(filename);

  size_t manager_i{0};
  for (const auto& factory_arg : dd.factory_args) {
    auto manager{make_structure_manager_stack_with_hypers_and_typeholder<
            ManagerTypeList_t>::apply(factory_arg["structure"],
                                      factory_arg["adaptors"])};
    std::cout << factory_arg["structure"]["filename"]<< std::endl;
    for (auto atom : manager) {
      std::cout << atom.get_atom_type()<<std::endl;
    }
    const auto & rep_infos{dd.ref_data.at("rep_info").template get<json>()};

    for (const auto & rep_info : rep_infos.at(manager_i)) {
      const auto & hypers = rep_info.at("hypers").template get<json>();
      const auto & ref_representation = rep_info.at("feature_matrix").template get<Std2DArray_t>();

      // for (auto& el : hypers.items()) {
      //   std::cout << el.key() << " : " << el.value() << "\n";
      // }
      RepresentationManagerSortedCoulomb<Manager_t> representation{manager,hypers};
      representation.compute();

      const auto & test_representation = representation.get_representation_full();

      for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {

        for (size_t col_i{0}; col_i < ref_representation[row_i].size();
              ++col_i) {
          auto diff{std::abs(ref_representation[row_i][col_i] -
                              test_representation(row_i, col_i))};
          if (diff > 1e-12) {
            std::cout << diff << "\n";
          }
        }
      }


    }


    manager_i += 1;
  }

  // using SMType =
  //     StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
  //                                AdaptorStrict>;

  // bool verbose{false};
  // // bool verbose_rep{false};
  // double cutoff{2.};
  // // bool consider_ghost_neighbours{false};
  // // std::string filename{"crystal_structure.json"};
  // std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};

  // // Initialize by hand
  // auto manager{make_structure_manager<StructureManagerCenters>()};
  // manager->update(filename);
  // auto pair_manager{
  //     make_adapted_manager<AdaptorNeighbourList>(manager, cutoff)};
  // pair_manager->update();
  // auto adaptor_strict{
  //     make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
  // adaptor_strict->update(filename);

  // // Or use the hyper thing
  // json adaptors_hypers = R"([
  //   {"name": "AdaptorNeighbourList", "initialization_arguments":{"cutoff": 2, "consider_ghost_neighbours": false}},
  //   {"name": "AdaptorStrict", "initialization_arguments":{"cutoff": 2}}
  // ])"_json;

  // json structure_inputs{{"filename", filename}};

  // json name_hypers = R"([
  //   {"name": "StrictNL", "initialization_arguments":{"cutoff": 2, "consider_ghost_neighbours": false}},
  //   {"name": "AdaptorStrict", "initialization_arguments":{"cutoff": 2}}
  // ])"_json;

  // auto my_man = wrap_factory<SMType>()(structure_inputs, adaptors_hypers);

  // std::cout << my_man->get_name() << std::endl;
  // auto lower_manager = extract_underlying_manager<-1>(my_man);
  // std::cout << lower_manager->get_name() << std::endl;

  // for (auto && center : my_man) {
  //   if (verbose) {
  //     std::cout << "################################# 2" << std::endl;
  //     std::cout << center.get_atom_type() << std::endl;
  //   }
  //   for (auto neigh : center) {
  //     if (verbose) {
  //       std::cout << neigh.get_atom_type() << ", ";
  //     }
  //   }
  //   if (verbose)
  //     std::cout << std::endl;
  // }

  // // Or use a fancier helper to do it in 1 line here
  // auto a1 = std::make_tuple(cutoff, false, cutoff);
  // auto a0 = std::make_tuple(filename);

  // using AdaptorTypeHolder_t = typename StructureManagerTypeHolder<
  //     StructureManagerCenters, AdaptorNeighbourList,
  //     AdaptorStrict>::type_list;
  // auto aa = std::make_tuple(a0, a1);
  // auto man{make_structure_manager_stack_with_tuple_and_typeholder<
  //     AdaptorTypeHolder_t>::apply(aa)};
  // std::cout << man->get_name() << std::endl;
  // // and there
  // auto mmmm = make_structure_manager_stack<StructureManagerCenters,
  //                                          AdaptorNeighbourList,
  //                                          AdaptorStrict>(
  //     filename, cutoff, consider_ghost_neighbours, cutoff);

  // MultipleStrictStructureManager<StructureManagerCenters> meta{};

  // for (auto && manager : meta.managers) {
  //   if (verbose) {
  //     std::cout << "################################# 1" << std::endl;
  //     std::cout << manager->size() << std::endl;
  //   }
  //   auto lower_manager = extract_underlying_manager<-2>(manager);
  //   std::cout << lower_manager->get_name() << std::endl;

  //   for (auto && center : manager) {
  //     if (verbose) {
  //       std::cout << center.get_atom_type() << std::endl;
  //       std::cout << "################################# 2" << std::endl;
  //     }
  //     for (auto neigh : center) {
  //       if (verbose) {
  //         std::cout << neigh.get_atom_type() << std::endl;
  //       }
  //     }
  //   }
  // }

  // json hypers{{"central_decay", 10},
  //               {"interaction_cutoff", 10},
  //               {"interaction_decay", 10},
  //               {"size", 20},
  //               {"sorting_algorithm", "distance"}};
  json hypers{{"interaction_cutoff", 6.0},
              {"cutoff_smooth_width", 1.0},
              {"max_radial", 10},
              {"max_angular", 8},
              {"gaussian_sigma_type", "Constant"},
              {"gaussian_sigma_constant", 0.5}};

  using Feature_t = FeatureManagerDense<double>;
  Feature_t features{810, hypers};
  // size_t i_center{0};
  // for (auto & manager : meta.managers) {
  //   // double central_decay{10};
  //   // double interaction_cutoff{10};
  //   // double interaction_decay{10};
  //   // size_t size{50};

  //   Representation_t representation{manager, hypers};
  //   representation.compute();
  //   features.insert(i_center, representation);
  //   i_center += manager->size();
  //   auto rep = representation.get_representation_full();
  //   if (verbose_rep) {
  //     std::cout << rep.size() << ", " << rep.cols() << ", " << rep.rows()
  //               << std::endl;
  //     for (auto ii{0}; ii < rep.cols(); ++ii) {
  //       for (auto jj{0}; jj < rep.rows(); ++jj) {
  //         std::cout << rep(jj, ii) << ", ";
  //       }
  //       std::cout << std::endl;
  //     }
  //   }
  // }

  // auto mat = features.get_feature_matrix();

  // for (size_t ii{0}; ii < mat.cols(); ii++) {
  //   for (size_t jj{0}; jj < mat.rows(); jj++) {
  //     std::cout << mat(jj, ii) << ", ";
  //   }
  //   std::cout << std::endl;
  // }

  return (0);
}
