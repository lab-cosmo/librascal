/**
 * file   test_nl.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2018
 *
 * @brief Example for Neighbour list
 *
 * Copyright © 2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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
#include "named_tuple.hh"
#include "enum_map.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

// using Representation_t = RepresentationManagerSortedCoulomb<
//     AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;
using Representation_t = RepresentationManagerSphericalExpansion<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

// template <class StructureManager>
// struct MultipleStrictStructureManager {
//   using Manager1_t = StructureManager;
//   using Manager2_t = AdaptorNeighbourList<Manager1_t>;
//   using Manager_t = AdaptorStrict<Manager2_t>;

//   MultipleStrictStructureManager() {
//     std::vector<std::string> filenames{{"alanine-X.json"}};
//     std::vector<double> cutoffs{{2, 3}};
//     bool consider_ghost_neighbours{false};
//     for (auto filename : filenames) {
//       for (auto cutoff : cutoffs) {
//         auto manager =
//             make_structure_manager_stack<StructureManager,
//             AdaptorNeighbourList,
//                                          AdaptorStrict>(
//                 filename, cutoff, consider_ghost_neighbours, cutoff);
//         this->managers.emplace_back(manager);
//       }
//     }
//   }

//   ~MultipleStrictStructureManager() {}

//   std::list<std::shared_ptr<Manager_t>> managers{};
// };




// template <size_t CurrentPosition, typename ManagerImplementation,
//           template <class> class AdaptorImplementation,
//           template <class> class... AdaptorImplementationPack>
// struct FactoryMap {
//   using Manager_t = AdaptorImplementation<ManagerImplementation>;
//   using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
//   using ManagerPtr_t = std::shared_ptr<Manager_t>;
//   constexpr static size_t NextPosition{CurrentPosition + 1};

//   using NextFactory_t = FactoryMap<NextPosition, Manager_t,
//                                               AdaptorImplementationPack...>;

//   //! General case
//   template <typename... Args, typename Hypers_t>
//   FactoryMap(ImplementationPtr_t & m,
//                         const Hypers_t & adaptors_hypers)
//       : manager{make_adapted_manager_hypers_util<
//             AdaptorImplementation,
//             CurrentPosition>::apply(m, adaptors_hypers)},
//         next_stack{manager, adaptors_hypers} {}

//   ManagerPtr_t manager;
//   NextFactory_t next_stack;

//   decltype(auto) get_manager() { return this->next_stack.get_manager(); }
// };

// //! End of recursion
// template <size_t CurrentPosition, typename ManagerImplementation,
//           template <class> class AdaptorImplementation>
// struct FactoryMap<CurrentPosition, ManagerImplementation,
//                               AdaptorImplementation> {
//   using Manager_t = AdaptorImplementation<ManagerImplementation>;
//   using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
//   using ManagerPtr_t = std::shared_ptr<Manager_t>;
//   using type = Manager_t;
//   constexpr static size_t NextPosition{CurrentPosition + 1};

//   template <typename... Args, typename Hypers_t>
//   FactoryMap(ImplementationPtr_t & m,
//                         const Hypers_t & adaptors_hypers)
//       : manager{make_adapted_manager_hypers_util<
//             AdaptorImplementation,
//             CurrentPosition>::apply(m, adaptors_hypers)} {}

//   ManagerPtr_t manager;

//   ManagerPtr_t get_manager() { return this->manager; }
// };


// template <typename ManagerImplementation,
//           template <class> class AdaptorImplementation,
//           template <class> class... AdaptorImplementationPack>
// struct AdaptorTypeList {
//   using Manager_t = AdaptorImplementation<ManagerImplementation>;
//   using type =
//       typename AdaptorTypeStacker<Manager_t,
//                                   AdaptorImplementationPack...>::type;
// };

// template <typename ManagerImplementation,
//           template <class> class AdaptorImplementation>
// struct AdaptorTypeList<ManagerImplementation, AdaptorImplementation> {
//   using type = AdaptorImplementation<ManagerImplementation>;
// };

template<typename StructureManagerTypeHolder>
decltype(auto) wrap_factory() {
  return []( const json& a, const json& b ) {
      return make_structure_manager_stack_with_hypers_and_typeholder<
                typename StructureManagerTypeHolder::type_list>::apply(a,b);
    };
  }

// decltype(auto) dispatch_factory(std::string name) {
//   using SMType1 = StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList>;
//   using SMType2 = StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList,AdaptorStrict>;
//   if (name == "one") {
//     return wrap_factory<SMType1>();
//   } else if (name == "two") {
//     return wrap_factory<SMType2>();
//   }
// }


enum class NeighbourList {LinkedCell, Strict};

constexpr std::initializer_list<
          std::pair<NeighbourList, const char*>> mapping = {
    {NeighbourList::LinkedCell, "LinkedCell"},
    {NeighbourList::Strict, "StrictNL"},
};

constexpr const char* get_name(NeighbourList e) {
    for (auto& p : mapping) {
        if (e == p.first) {
            return p.second;
        }
    }
}
constexpr NeighbourList get_e(const char* name) {
    for (auto& p : mapping) {
        if (name == p.second) {
            return p.first;
        }
    }
}

template<NeighbourList NL>
decltype(auto) dispatch_factory();

template <>
decltype(auto) dispatch_factory<NeighbourList::Strict>() {
  using SMType2 = StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList,AdaptorStrict>;
  return wrap_factory<SMType2>();
}

int main() {

  std::cout << get_name(NeighbourList::LinkedCell) << std::endl;
  std::cout << static_cast<int>(get_e("StrictNL")) << std::endl;

  using TypeList = StructureManagerTypeList<
    StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList>,
    StructureManagerTypeHolder<
            StructureManagerCenters,AdaptorNeighbourList,AdaptorStrict>
              >;
  using SMType1 = StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList>;
  using SMType2 = StructureManagerTypeHolder<StructureManagerCenters,AdaptorNeighbourList,AdaptorStrict>;

  auto name{
      internal::GetTypeNameHelper<typename SMType1::type>::GetTypeName()};
  std::cout << name;
  // std::cout << const_hash<name>();
  auto factory_map1 = make_named_tuple(
      make_named_value(param("LinkedCell"), wrap_factory<SMType1>())
    );

  auto factory_map2 = make_named_tuple(
      make_named_value(param("StrictNL"), wrap_factory<SMType2>())
    );

  auto factory_map = named_tuple_cat(factory_map1, factory_map2);

  bool verbose{false};
  bool verbose_rep{false};
  double cutoff{2.};
  bool consider_ghost_neighbours{false};
  // std::string filename{"crystal_structure.json"};
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};

  // Initialize by hand
  auto manager{make_structure_manager<StructureManagerCenters>()};
  manager->update(filename);
  auto pair_manager{
      make_adapted_manager<AdaptorNeighbourList>(manager, cutoff)};
  pair_manager->update();
  auto adaptor_strict{
      make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
  adaptor_strict->update(filename);

  // Or use the hyper thing
  json adaptors_hypers = R"([
    {"name": "AdaptorNeighbourList", "initialization_arguments":{"cutoff": 2, "consider_ghost_neighbours": false}},
    {"name": "AdaptorStrict", "initialization_arguments":{"cutoff": 2}}
  ])"_json;


  json structure_inputs{{"filename", filename}};
  // constexpr const char* name_ = "StrictNL";
  // auto my_man = factory_map[param(name_)](structure_inputs, adaptors_hypers);
  json name_hypers = R"([
    {"name": "StrictNL", "initialization_arguments":{"cutoff": 2, "consider_ghost_neighbours": false}},
    {"name": "AdaptorStrict", "initialization_arguments":{"cutoff": 2}}
  ])"_json;

  const auto E{get_name(aa)};
  auto my_man = dispatch_factory<E>()(structure_inputs, adaptors_hypers);

  std::cout << my_man->get_name() << std::endl;
  auto lower_manager = extract_underlying_manager<-1>(my_man);
  std::cout << lower_manager->get_name() << std::endl;

  for (auto && center : my_man) {
    if (verbose) {
      std::cout << "################################# 2" << std::endl;
      std::cout << center.get_atom_type() << std::endl;
    }
    for (auto neigh : center) {
      if (verbose) {
        std::cout << neigh.get_atom_type() << ", ";
      }
    }
    if (verbose) std::cout << std::endl;
  }

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
  size_t i_center{0};
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
