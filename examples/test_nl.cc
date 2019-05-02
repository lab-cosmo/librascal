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


struct TestData {
  using ManagerTypeHolder_t =
      StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
                                 AdaptorStrict>;

  TestData() = default;

  void get_ref(const std::string & ref_filename) {
    std::vector<std::uint8_t> ref_data_ubjson;
    internal::read_binary_file(ref_filename, ref_data_ubjson);
    this->ref_data = json::from_ubjson(ref_data_ubjson);
    auto filenames =
        this->ref_data.at("filenames").get<std::vector<std::string>>();
    auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

    for (auto && filename : filenames) {
      for (auto && cutoff : cutoffs) {
        json parameters;
        json structure{{"filename", filename}};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
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
      StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
                                 AdaptorStrict>;
  using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
  using Manager_t = typename ManagerTypeHolder_t::type;
  using Std2DArray_t = std::vector<std::vector<double>>;

  auto dd{TestData()};
  std::string filename{"reference_data/sorted_coulomb_reference.ubjson"};
  dd.get_ref(filename);

  size_t manager_i{0};
  for (const auto & factory_arg : dd.factory_args) {
    auto manager{make_structure_manager_stack_with_hypers_and_typeholder<
        ManagerTypeList_t>::apply(factory_arg["structure"],
                                  factory_arg["adaptors"])};
    std::cout << factory_arg["structure"]["filename"] << std::endl;
    for (auto atom : manager) {
      std::cout << atom.get_atom_type() << std::endl;
    }
    const auto & rep_infos{dd.ref_data.at("rep_info").template get<json>()};

    for (const auto & rep_info : rep_infos.at(manager_i)) {
      const auto & hypers = rep_info.at("hypers").template get<json>();
      const auto & ref_representation =
          rep_info.at("feature_matrix").template get<Std2DArray_t>();

      RepresentationManagerSortedCoulomb<Manager_t> representation{manager,
                                                                   hypers};
      representation.compute();

      const auto & test_representation =
          representation.get_representation_full();

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


  return (0);
}
