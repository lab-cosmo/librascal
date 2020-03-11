/**
 * @file   performance/profiles/profile_sparse_kernel.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   12 March 2020
 *
 * @brief  Example for profiling the spherical invariants (SOAP)
 *
 * Copyright © 2020 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/models/sparse_kernels.hh"
#include "rascal/models/pseudo_points.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils/utils.hh"

#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <string>

// using namespace std;
using namespace rascal;  // NOLINT

using ManagerTypeHolder_t =
    StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
                               AdaptorCenterContribution, AdaptorStrict>;
using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
using Manager_t = typename ManagerTypeHolder_t::type;
using ManagerCollection_t =
    typename TypeHolderInjector<ManagerCollection, ManagerTypeList_t>::type;
using Representation_t = CalculatorSphericalInvariants;

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide atomic structure json filename as argument";
    std::cerr << std::endl;
    return -1;
  }

  json data = json_io::load(argv[1]);

  const int N_ITERATIONS = data["N_ITERATIONS"].get<int>();
  const int N_SPARSE_POINTS = data["N_SPARSE_POINTS"].get<int>();
  json hypers = data["rep_hypers"].get<json>();
  json kernel_hypers = data["kernel_hypers"].get<json>();

  std::string filename{data["filename"].get<std::string>()};
  double cutoff{data["cutoff"].get<double>()};
  int begin{data["start"].get<int>()};
  int length{data["length"].get<int>()};
  json structure{{"filename", filename}};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments", {{"cutoff", cutoff}, {"skin", 0.}}}};
  json ad1b{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad1b);
  adaptors.emplace_back(ad2);

  ManagerCollection_t collection{adaptors};
  ManagerCollection_t managers{adaptors};
  Representation_t soap{hypers};
  SparseKernel kernel{kernel_hypers};

  collection.add_structures(filename, begin, length);
  soap.compute(collection);
  auto man = collection[0];
  managers.add_structure(man);

  PseudoPointsBlockSparse<Representation_t> pseudo_points{};

  std::vector<std::vector<int>> selected_ids;
  int n_centers{0};
  for (auto & manager : collection) {
    selected_ids.emplace_back();
    int ii{0};
    for (auto center : manager) {
      if (n_centers >= N_SPARSE_POINTS) {
        break;
      }
      selected_ids.back().push_back(ii);
      ++ii;
      ++n_centers;
    }
  }
  pseudo_points.push_back(soap, collection, selected_ids);

  // auto KNM_test{kernel.compute(representation, managers, pseudo_points)};

  std::cout << "structure filename: " << filename << std::endl;

  std::chrono::duration<double> elapsed{};

  auto start = std::chrono::high_resolution_clock::now();
  // This is the part that should get profiled
  for (int looper{0}; looper < N_ITERATIONS; looper++) {
    auto KNM_der{kernel.compute_derivative(soap, managers, pseudo_points)};
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "gradient kernel"
            << " elapsed: " << elapsed.count() / (N_ITERATIONS* man->size())
            << " seconds / atom" << std::endl;
}
