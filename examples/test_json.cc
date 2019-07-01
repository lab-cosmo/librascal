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
#include "json_io.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

using Representation_t = RepresentationManagerSOAP<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;


using LayerByOrder = std::index_sequence<0, 0, 1>;





// template<std::size_t... I>
// constexpr auto sum_seq(std::index_sequence<I...>, const int& sum = 0) {
//   return std::make_tuple(a[I]...);
// }

// template<size_t I, size_t... Is,size_t Sum = 0 >
// constexpr auto sum(const int& sum) {
//   return sum + 1 + I;
// }

int main() {
  double cutoff{3.};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments",
            {{"cutoff", cutoff}, {"consider_ghost_neighbours", false}}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);

  json j1{{1,23,4,6,7},{1,23,4,6,9},{4,23,4,6,9}};
  std::cout << j1.dump() <<std::endl;

  auto mat1 = j1.get<Eigen::MatrixXd>();

  std::cout << mat1 <<std::endl;

  AtomicStructure<3> sss{};
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  sss.set_structure(filename);

  std::cout << sss.atom_types.rows() << ", "<<sss.atom_types.cols() << std::endl;
  std::cout << sss.cell.rows() << ", "<<sss.cell.cols() << std::endl;
  std::cout << sss.positions.rows() << ", "<<sss.positions.cols() << std::endl;
  std::cout << sss.pbc.rows() << ", "<<sss.pbc.cols() << std::endl;

  std::cout << sss.cell << std::endl;
  std::cout << sss.positions << std::endl;

  auto manager = make_structure_manager<StructureManagerCenters>();
  manager->update(sss);
  // auto mat2 = j1.get<std::vector<double>>();
  // for (auto& el : mat2) {
  //   std::cout << el;
  // }
  // std::cout <<std::endl;
  std::cout << internal::SumSequence<LayerByOrder>::value << std::endl;

  auto man = stack_adaptors<StructureManagerCenters, AdaptorNeighbourList, AdaptorStrict>(manager, adaptors);

  std::cout << man->get_name() << std::endl;
  std::cout << extract_underlying_manager<2>(man)->get_name() << std::endl;
  std::cout << extract_underlying_manager<1>(man)->get_name() << std::endl;
  auto man1 = extract_underlying_manager<0>(man);
  std::cout << man1->get_name() << std::endl;

  return (0);
}
