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
 * Copyright  2018 markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include <structure_managers/adaptor_half_neighbour_list.hh>
#include <structure_managers/adaptor_increase_maxorder.hh>
#include <structure_managers/adaptor_neighbour_list.hh>
#include <structure_managers/adaptor_strict.hh>
#include <structure_managers/make_structure_manager.hh>
#include <structure_managers/property.hh>
#include <structure_managers/structure_manager_centers.hh>

/**
 * This small example highlight the possibilities of adapting a structure
 * (e.g. read in from a file) to certain iterations. By stacking `Manager` and
 * `Adaptors`, the necessary strictness of a neighbour list or maximum Order
 * (e.g. triplets) can be constructed. It is all backwards compatible, i.e. one
 * can always iterate over atoms, but then increasingly higher Orders, depending
 * on the stacking.
 */

// Shorthands for templated types for readability.
using Manager_t = rascal::StructureManagerCenters;
using PairManager_t = rascal::AdaptorNeighbourList<Manager_t>;

using StrictPairManager_t = rascal::AdaptorStrict<PairManager_t>;
using TripletManager_t = rascal::AdaptorMaxOrder<StrictPairManager_t>;

int main() {
  auto manager{rascal::make_structure_manager<Manager_t>()};
  double cutoff{1.};
  /**
   * These structures here are sample structures and can be used to iterate
   * over.
   *
   * `crystal_structure.json` is a fully periodic metallic hcp structure.
   *
   * `alanine-X.json` is a unit cell of a polyalanine chain.
   *
   * `simple_cubic_9.json` is an artificial 9-atom test structure.
   */

  std::string filename{"crystal_structure.json"};
  // std::string filename{"alanine-X.json"};
  // std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};

  std::cout << "Reading structure " << filename << std::endl;

  // `manager` is the data object which reads, stores and gives access to atom
  // positions, types. It also provides iteration over all atom.
  //  manager.update(filename);

  // `pair_manager` is constructed with the `manager` and a `cutoff`.
  auto pair_manager{rascal::make_adapted_manager<rascal::AdaptorNeighbourList>(
      manager, cutoff, true)};
  // By invoking the `.update()` method, a neighbour list is built.
  //  pair_manager->update();

  // `strict_manager` is constructed with a `pair_manager`.
  auto strict_manager{rascal::make_adapted_manager<rascal::AdaptorStrict>(
      pair_manager, cutoff)};
  // calling the `.update()` method triggers the build of a strict neighbourlist
  // (all pairs are within the specified cutoff)
  //  strict_manager.update();

  // `triplet_manager` is constructed with a pair list (strict or not, here
  // strict)
  auto triplet_manager{
      rascal::make_adapted_manager<rascal::AdaptorMaxOrder>(strict_manager)};
  // `.update()` triggers the extension of the pair list to triplets
  // triplet_manager->update(positions, atom_types, cell, PBC_t{pbc.data()});
  triplet_manager->update(filename);

  // Iteration over `manager`
  std::cout << "manager iteration over atoms" << std::endl;
  for (auto atom : manager) {
    std::cout << "atom " << atom.get_atom_tag() << " global index "
              << atom.get_global_index() << std::endl;
  }

  // `pair_manager` provides iteration over atoms and pairs
  for (auto atom : pair_manager) {
    for (auto pair : atom) {
      std::cout << "pair (" << atom.get_atom_tag() << ", "
                << pair.get_atom_tag() << " ) global index "
                << pair.get_global_index() << std::endl;
    }
  }

  // `strict_manager` provides iteration over atoms and strict pairs
  for (auto atom : strict_manager) {
    for (auto pair : atom) {
      std::cout << "strict pair (" << atom.get_atom_tag() << ", "
                << pair.get_atom_tag() << ") global index "
                << pair.get_global_index() << std::endl;
    }
  }

  // `triplet_manager` provides iteration over atoms, strict pairs and strict
  // triplets
  for (auto atom : triplet_manager) {
    for (auto pair : atom) {
      for (auto triplet : pair) {
        std::cout << "triplet (" << atom.get_atom_tag() << ", "
                  << pair.get_atom_tag() << ", " << triplet.get_atom_tag()
                  << ") global index " << triplet.get_global_index()
                  << std::endl;
      }
    }
  }
}
