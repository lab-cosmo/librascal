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

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>

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

template<char x, char... xs>
struct hash_calc {
    static constexpr std::uint64_t apply () {
       return  (hash_calc<xs...>::apply() ^ x) * 16777619u;
    };
};

template<char x>
struct hash_calc<x> {
    static constexpr std::uint64_t apply () {
       return  2166136261u;
    };
};

template<char... xs>
constexpr std::uint64_t hash () {
    return hash_calc<xs...>::apply();
}

// only clang/gcc compatible
template <typename Char, Char... Cs>
constexpr auto operator""_c() {
  // constexpr auto hash_{foonathan::string_id::detail::sid_hash(s)};
  // return std::integral_constant<std::size_t, hash<Cs...>()>{};
  return fn_detail::make_named_param< std::integral_constant<foonathan::string_id::detail::hash_type, hash<Cs...>()> >{};
};


template <typename Tuple1, size_t... Indices1, typename Tuple2, size_t... Indices2>
decltype(auto) tuple_cat1(Tuple1&& tup1, Tuple2&& tup2,
                std::index_sequence<Indices1...>, std::index_sequence<Indices2...>)
{
  auto aa = make_named_tuple(
    get<Indices1>(std::forward<Tuple1>(tup1))...,
    get<Indices2>(std::forward<Tuple2>(tup2))...
  );
  // std::cout << aa;
  return aa;
}

template< class T >
class named_tuple_size;

template< class... Types >
class named_tuple_size< fn_detail::named_tuple<Types...> >
  : public std::integral_constant<std::size_t, sizeof...(Types)> { };


template <typename Tuple1, typename Tuple2>
decltype(auto) named_tuple_cat(Tuple1&& tup1, Tuple2&& tup2)
{
  tuple_cat1(
   std::forward<Tuple1>(tup1),
   std::forward<Tuple2>(tup2),
   std::make_index_sequence<named_tuple_size<std::decay_t<Tuple1>>::value>{},
   std::make_index_sequence<named_tuple_size<std::decay_t<Tuple2>>::value>{}
  );
}



int main() {


  // auto tup1 = make_named_tuple("first"_c = 3, "second"_c = 12.7);
  // auto first = tup1["first"_c];
  // auto second = tup1["second"_c];
  // std::cout << first << ", " << second << std::endl;
  // auto tup2 = make_named_tuple("third"_c = 3.14);
  // std::cout << tup2.get_by_index<0>()  << std::endl;

  // // auto tup3 = named_tuple_cat(tup1, tup2);
  // // auto third = tup3["third"_c];
  // // // fn_detail::named_tuple<fn_detail::named_param<std::integral_constant<long unsigned int, 10726708487119078247>, int> , fn_detail::named_param<
  // std::cout << third  << std::endl;

  auto func1 = []( auto a, auto b ) {
          return make_structure_manager_stack<
          StructureManagerCenters,AdaptorNeighbourList>(a,b);
        };

  auto factory_map1 = make_named_tuple(
      // "first"_c = aa,
      "first"_c = func1
    );

  // std::cout << factory_map1;
  auto func = []( auto a, auto b ) {
          return make_structure_manager_stack<
          StructureManagerCenters,AdaptorNeighbourList,AdaptorStrict>(a,b);
        };

  auto factory_map2 = make_named_tuple(
      // "first"_c = aa,
      "second"_c = func
    );

  // auto factory_map = named_tuple_cat(factory_map1, factory_map2);
  // std::cout << named_tuple_cat(factory_map1, factory_map2);

  bool verbose{true};
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

  // auto my_man =
  //     make_structure_manager_stack<StructureManagerCenters,
  //                                  AdaptorNeighbourList, AdaptorStrict>(
  //         structure_inputs, adaptors_hypers);
  // const auto second = factory_map2["second"_c];
  auto my_man = factory_map2["second"_c](structure_inputs, adaptors_hypers);
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
