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

#include "basic_types.hh"
#include "models/kernels.hh"
#include "rascal_utility.hh"
#include "representations/calculator_sorted_coulomb.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"
#include "structure_managers/adaptor_center_contribution.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/make_structure_manager.hh"
#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/structure_manager_collection.hh"

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <random>
#include <string>

using namespace rascal;  // NOLINT

int main() {
  // Test1()();
  // std::string filename{"reference_data/dft-smiles_500.ubjson"};
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  //  std::string filename{"reference_data/alloy-small.json"};
  // std::string filename{"reference_data/alloy-small.json"};
  // std::string filename{"reference_data/diamond_cubic.json"};
  std::string rep_id{"pp"};

  double cutoff{3.};

  json structure{{"filename", filename}};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments",
            {{"cutoff", cutoff},
             {"consider_ghost_neighbours", false},
             {"skin", 0.}}}};

  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);
  /*auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structure, adaptors);
    */
  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>(
          structure, adaptors);

  std::cout << "n_centers: " << manager->size() << std::endl;
  for (auto center : manager) {
    auto ctag = center.get_atom_tag();
    std::cout << "Center: " << ctag << " n. neighbors " << center.size()
              << std::endl;

    for (auto neigh : center) {
      auto tag_list = neigh.get_atom_tag_list();

      auto atom_j = neigh.get_atom_j();
      auto atom_j_tag = atom_j.get_atom_tag_list();
      auto atom_j_ids = atom_j.get_cluster_indices();
      std::cout << "neigh: " << tag_list[0] << ", " << tag_list[1] << ", "
                << " tag_j: " << atom_j_tag[0] << ", " << atom_j_ids[0]
                << std::endl;
    }
  }

  return (0);
}
