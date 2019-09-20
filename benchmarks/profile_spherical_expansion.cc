/**
 * file   spherical_expansion_profile.cc
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   7 May 2019
 *
 * @brief  Example for profiling the spherical expansion
 *
 * Copyright © 2018 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/make_structure_manager.hh"
#include "rascal_utility.hh"
#include "representations/calculator_sorted_coulomb.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"
#include "basic_types.hh"
#include "atomic_structure.hh"

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>
#include <chrono>

// using namespace std;
using namespace rascal;  // NOLINT

const int N_ITERATIONS = 1000;

using Representation_t = CalculatorSphericalInvariants;
using Manager_t = AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>;
using Prop_t = typename CalculatorSphericalInvariants::Property_t<Manager_t>;
using PropDer_t =
    typename CalculatorSphericalInvariants::PropertyGradient_t<Manager_t>;

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide atomic structure json filename as argument";
    std::cerr << std::endl;
    return -1;
  }

  // TODO(max) put these in a file so they can be varied systematically
  // maybe together with the filename and iteration count
  std::string filename{argv[1]};

  double cutoff{5.};
  json hypers{
      {"max_radial", 8}, {"max_angular", 6}, {"compute_gradients", false}};

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
            {{"cutoff", cutoff}, {"consider_ghost_neighbours", false}}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);
  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>(
          structure, adaptors);

  AtomicStructure<3> ast{};
  ast.set_structure(filename);

  std::cout << "structure filename: " << filename << std::endl;

  std::chrono::duration<double> elapsed{};

  auto start = std::chrono::high_resolution_clock::now();
  // This is the part that should get profiled
  for (size_t looper{0}; looper < N_ITERATIONS; looper++) {
    manager->update(ast);
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Neighbour List"
            << " elapsed: " << elapsed.count() / N_ITERATIONS << " seconds"
            << std::endl;

  Representation_t representation{hypers};

  start = std::chrono::high_resolution_clock::now();
  // This is the part that should get profiled
  for (size_t looper{0}; looper < N_ITERATIONS; looper++) {
    representation.compute(manager);
  }
  finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Compute representation"
            << " elapsed: " << elapsed.count() / N_ITERATIONS << " seconds"
            << std::endl;

  // auto expn = representation.get_representation_full();
  // std::cout << "Sample SphericalExpansion elements " << std::endl
  //           << expn(0, 0) << " " << expn(0, 1) << " " << expn(0, 2) << "\n"
  //           << expn(1, 0) << " " << expn(1, 1) << " " << expn(1, 2) << "\n"
  //           << expn(2, 0) << " " << expn(2, 1) << " " << expn(2, 2) << "\n";

  // Profile again, this time with gradients
  hypers["compute_gradients"] = true;
  Representation_t representation_gradients{hypers};
  start = std::chrono::high_resolution_clock::now();
  // This is the part that should get profiled
  for (size_t looper{0}; looper < N_ITERATIONS; looper++) {
    representation_gradients.compute(manager);
  }
  finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_grad{};
  elapsed_grad = finish - start;
  std::cout << "Compute representation with gradients"
            << " elapsed: " << elapsed_grad.count() / N_ITERATIONS << " seconds"
            << std::endl;
  std::cout << "Ratio (with gradients / without gradients): "
            << elapsed_grad.count() / elapsed.count() << std::endl;

  // auto expn2 = representation_gradients.get_representation_full();
  // std::cout << "Sample SphericalExpansion elements (should be identical) "
  //           << std::endl
  //           << expn2(0, 0) << " " << expn2(0, 1) << " " << expn2(0, 2) <<
  //           "\n"
  //           << expn2(1, 0) << " " << expn2(1, 1) << " " << expn2(1, 2) <<
  //           "\n"
  //           << expn2(2, 0) << " " << expn2(2, 1) << " " << expn2(2, 2) <<
  //           "\n";
  // TODO(max) print out analogous gradient components, for now see
  // spherical_expansion_example
}
