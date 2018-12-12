/**
 * file   structure_manager_iteration.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   22 Nov 2018
 *
 * @brief highlights the building and iteration possibilities of the concept
 *        'structure manager' and adaptors
 *
 * Copyright © 2018 markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include <structure_managers/structure_manager_centers.hh>
#include <structure_managers/adaptor_neighbour_list.hh>
#include <structure_managers/adaptor_half_neighbour_list.hh>
#include <structure_managers/adaptor_strict.hh>
#include <structure_managers/adaptor_increase_maxorder.hh>
#include <structure_managers/property.hh>

/**
 * Shorthands for templated types.
 */
using Manager_t = rascal::StructureManagerCenters;
using PairManager_t = rascal::AdaptorNeighbourList<Manager_t>;
using PairManagerHalf_t = rascal::AdaptorHalfList<PairManager_t>;

using StrictPairManager_t = rascal::AdaptorStrict<PairManager_t>;
using TripletManager_t = rascal::AdaptorMaxOrder<StrictPairManager_t>;
using TripletManager2_t = rascal::AdaptorMaxOrder<PairManager_t>;

int main() {
  Manager_t manager;
  double cutoff{1.};
  // std::string filename{"crystal_structure.json"};
  // std::string filename{"alanine-X.json"};
  std::string filename{"simple_cubic_9.json"};

  std::cout << "filename " << filename << std::endl;

  manager.update(filename);

  PairManager_t pair_manager{manager, cutoff};
  pair_manager.update();

  for (auto atom : pair_manager) {
    for (auto pair : atom) {
      std::cout << "pair (" << atom.get_atom_index() << ", "
                << pair.get_atom_index() << " )" << std::endl;
    }
  }

  // PairManagerHalf_t pair_manager_half{pair_manager};
  // pair_manager_half.update();

  std::cout << "Building strict manager" << std::endl;
  StrictPairManager_t strict_manager{pair_manager, cutoff};
  strict_manager.update();

  for (auto atom : strict_manager) {
    for (auto pair : atom) {
      std::cout << "strict pair " << atom.get_atom_index() << ", "
                << pair.get_atom_index() << " global index "
                << pair.get_global_index() << std::endl;
    }
  }

  std::cout << "Building triplet manager" << std::endl;
  TripletManager_t triplet_manager{strict_manager};
  triplet_manager.update();

  std::cout << "Iteration over triplet manager " << std::endl;

  for (auto atom : triplet_manager) {
    std::cout << "atom " << atom.get_atom_index() << std::endl;

    for (auto pair : atom) {
      std::cout << "pair " << pair.get_atom_index() << std::endl;

      for (auto triplet : pair) {
        std::cout << "triplet (" << atom.get_atom_index() << ", "
                  << pair.get_atom_index() << ", " << triplet.get_atom_index()
                  << ")" << std::endl;
      }
    }
  }
}
