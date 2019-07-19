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

// #include "structure_managers/structure_manager_centers.hh"
// #include "structure_managers/adaptor_strict.hh"
// #include "structure_managers/adaptor_neighbour_list.hh"
// #include "structure_managers/make_structure_manager.hh"
// #include "rascal_utility.hh"
// #include "representations/calculator_sorted_coulomb.hh"
// #include "representations/calculator_spherical_expansion.hh"
// #include "representations/calculator_spherical_invariants.hh"
// #include "representations/feature_manager_dense.hh"
// #include "representations/feature_manager_block_sparse.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <vector>
#include <memory>
#include <initializer_list>
// #include <typeinfo>
// #include <cstdlib>
// #include <cxxabi.h>

// using namespace std;
using namespace rascal;  // NOLINT

// std::string demangle(const char* name) {

//   int status = -4; // some arbitrary value to eliminate the compiler warning

//   // enable c++11 by passing the flag -std=c++11 to g++
//   std::unique_ptr<char, void(*)(void*)> res {
//       abi::__cxa_demangle(name, NULL, NULL, &status),
//       std::free
//   };

//   return (status==0) ? res.get() : name ;
// }

// struct Base {
//   virtual ~Base() = default;
// };

// template<size_t Num>
// struct Impl : Base {
//   static const size_t value = Num;
//   const std::string STRING = "other";
// };

// template<>
// struct Impl<1> : Base {
//   static const size_t value = 1;
//   const std::string STRING = "one";
// };

// template<>
// struct Impl<2> : Base {
//   static const size_t value = 2;
//   const std::string STRING = "two";
// };

// template <size_t Num>
// decltype(auto) make() {
//   return std::static_pointer_cast<Base>(
//       std::make_shared<Impl<Num>>());
// }

// template <size_t Num>
// decltype(auto) downcast(
//     const std::shared_ptr<Base> & base) {
//   return std::dynamic_pointer_cast<Impl<Num>>(base);
// }

// template <class Prop, class Manager>
// decltype(auto) make_prop(Manager& manager) {
//   return std::static_pointer_cast<PropertyBase>(
//       std::make_shared<Prop>(*manager));
// }

// template <class Prop>
// decltype(auto) downcast_prop(
//     const std::shared_ptr<PropertyBase> & base) {
//   return std::dynamic_pointer_cast<Prop>(base);
// }

// using Manager_t =
// AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>; using
// Property1_t = BlockSparseProperty<double, 1, 0, Manager_t, std::vector<int>>;
// using Property2_t = BlockSparseProperty<double, 1, 0, Manager_t,
// std::vector<double>>;
int main() {

  // std::vector<std::shared_ptr<Base>> tests{};
  // tests.emplace_back(make<1>());
  // tests.emplace_back(make<3>());
  // tests.emplace_back(make<2>());
  // tests.emplace_back(make<10>());

  // auto res0{downcast<1>(tests[0])};
  // std::cout << res0->STRING << std::endl;

  // auto res1{downcast<3>(tests[1])};
  // std::cout << res1->STRING << std::endl;
  // // std::cout << demangle(typeid(res1).name()) << std::endl;

  // double cutoff{3.};
  // std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  // json structure{{"filename", filename}};
  // json adaptors;
  // json ad1{{"name", "AdaptorNeighbourList"},
  //          {"initialization_arguments",
  //                   {{"cutoff", cutoff}, {"consider_ghost_neighbours",
  //                   false}}}};
  // json ad2{{"name", "AdaptorStrict"},
  //          {"initialization_arguments", {{"cutoff", cutoff}}}};
  // adaptors.emplace_back(ad1);
  // adaptors.emplace_back(ad2);
  // auto manager =
  //         make_structure_manager_stack<StructureManagerCenters,
  //                 AdaptorNeighbourList, AdaptorStrict>(
  //                 structure, adaptors);

  // auto prop1_base{make_prop<Property1_t>(manager)};
  // std::cout << internal::GetTypeName<Property1_t>() << std::endl;
  // std::cout << internal::GetTypeNameHelper<Property1_t>::GetTypeName() <<
  // std::endl; std::cout << internal::GetTypeNameHelper<double>::GetTypeName()
  // << std::endl; Property1_t::check_compatibility(*prop1_base);

  // Property2_t::check_compatibility(*prop1_base);

  // auto prop1{downcast_prop<Property1_t>(prop1_base)};

  return (0);
}
