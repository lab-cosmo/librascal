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
#include "basic_types.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>

// using namespace std;
using namespace rascal;  // NOLINT

// using Manager_t = StructureManagerCenters;
constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

using Representation_t = RepresentationManagerSortedCoulomb<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

template <class StructureManager>
struct MultipleStrictStructureManager {
  using Manager1_t = StructureManager;
  using Manager2_t = AdaptorNeighbourList<Manager1_t>;
  using Manager_t = AdaptorStrict<Manager2_t>;

  MultipleStrictStructureManager() {
    std::vector<std::string> filenames{
      {"alanine-X.json"}
                                        };
    std::vector<double> cutoffs{{3, 4}};
    bool consider_ghost_neighbours{false};
    for (auto filename : filenames) {
      for (auto cutoff : cutoffs) {
        // auto manager = make_structure_manager_stack<
        //       StructureManager,AdaptorNeighbourList,AdaptorStrict>
        //         (filename, std::make_tuple(cutoff), std::make_tuple(cutoff));
        auto manager = make_structure_manager_stack<
              StructureManager,AdaptorNeighbourList,AdaptorStrict>
                (filename, cutoff, consider_ghost_neighbours, cutoff);
        this->managers.emplace_back(
                manager
                );
      }
    }
  }

  ~MultipleStrictStructureManager() {}

  // std::list<std::shared_ptr<Manager1_t>> managers1{};
  // std::list<std::shared_ptr<Manager2_t>> managers2{};
  std::list<std::shared_ptr<Manager_t>> managers{};
};


// template<template<class...>class U, typename MI, typename typeholder>
// struct test;

// template<template<class...>class U, typename MI, template<class> class ...T>
// struct test<U, MI, AdaptorTypeHolder<T...>> {
//   template<typename ...Args>
//   static decltype(auto) func(Args ...args) {
//     return U<MI, T...>::apply(args...);
//   }
// };



// template<typename TemplateTypeHolder>
// struct call_with_typeholders;

// template<typename ...T>
// struct call_with_typeholders<std::tuple<T...>> {

//   template<typename ...Args>
//   static decltype(auto) make_manager_stack(Args ...args) {
//     return TypeExtractor<T...>::apply(args...);
//   }

//  protected:


// };

// template<typename TemplateTypeHolder, typename ArgsTypeHolder>
// struct call_with_typeholders;

// template<typename ...T, typename ...Args>
// struct call_with_typeholders<std::tuple<T...>, std::tuple<Args...>> {

//   static decltype(auto) make_manager_stack(std::tuple<Args...> tuple) {
//     return TypeExtractor<T...>::apply(tuple);
//   }

//  protected:

//   template<typename MI, typename TemplateTypeHolder_>
//   struct TypeExtractor;

//   template<typename MI, template<class> class ...Ti>
//   struct TypeExtractor<MI, AdaptorTypeHolder<Ti...>> {
//     using Manager_t = typename internal::AdaptorTypeStacker<MI,Ti...>::type;
//     using ManagerPtr_t = std::shared_ptr<Manager_t>;


//   };
// };


int main() {
  bool verbose{false};
  bool verbose_rep{false};
  double cutoff{2.};
  // std::string filename{"crystal_structure.json"};
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};

  auto manager{make_structure_manager<StructureManagerCenters>()};
  manager->update(filename);
  auto pair_manager{make_adapted_manager<AdaptorNeighbourList>(manager, cutoff)};
  pair_manager->update();
  auto adaptor_strict{make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
  adaptor_strict->update(filename);






  auto a1 = std::make_tuple(cutoff,false, cutoff);
  auto a0 = std::make_tuple(filename);
  // using AdaptorTypeHolder_t = AdaptorTypeHolder<AdaptorNeighbourList, AdaptorStrict>;
  //using Factory_t = std::tuple<std::string,std::tuple<double>,std::tuple<double>>;
  using AdaptorTypeHolder_t = typename StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList, AdaptorStrict>::type_list;
  auto aa = std::make_tuple(a0, a1);
  auto man{make_structure_manager_stack_with_tuple_and_typeholder<AdaptorTypeHolder_t>::apply(aa)};
  std::cout << man->get_name()<< std::endl;


  MultipleStrictStructureManager<StructureManagerCenters> meta{};

  for (auto&& manager : meta.managers) {
    // manager->update("alanine-X.json");
    // manager->update(positions, atom_types, cell, PBC_t{pbc.data()});
    if (verbose) {
      std::cout << "################################# 1"<< std::endl;
      std::cout << manager->size() << std::endl;
    }
    auto lower_manager = extract_underlying_manager<-2>(manager);
    std::cout << lower_manager->get_name()<< std::endl;

    for (auto&& center : manager) {
      if (verbose) {
        std::cout << center.get_atom_type() << std::endl;
        std::cout << "################################# 2"<< std::endl;
      }
      for (auto neigh : center) {
        if (verbose) {
          std::cout << neigh.get_atom_type() << std::endl;
        }
      }
    }
  }

  for (auto & manager : meta.managers) {
    // double central_decay{10};
    // double interaction_cutoff{10};
    // double interaction_decay{10};
    // size_t size{50};
    json hypers{{"central_decay", 10},
                {"interaction_cutoff", 10},
                {"interaction_decay", 10},
                {"size", 50},
                {"sorting_algorithm", "distance"}};
    Representation_t representation{manager, hypers};
    representation.compute();

    auto rep = representation.get_representation_full();
    if (verbose_rep) {
      std::cout << rep.size() << ", " << rep.cols() << ", " << rep.rows()
                << std::endl;
      for (auto ii{0}; ii < rep.cols(); ++ii) {
        for (auto jj{0}; jj < rep.rows(); ++jj) {
          std::cout << rep(jj, ii) << ", ";
        }
        std::cout << std::endl;
      }
    }
  }

  return(0);
}
