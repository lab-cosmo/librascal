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
#include <structure_managers/make_structure_manager.hh>
#include <structure_managers/adaptor_increase_maxorder.hh>
#include <structure_managers/property.hh>

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
  Eigen::MatrixXd positions(22, 3);
  Eigen::VectorXi atom_types(22);
  Eigen::MatrixXd cell(3, 3);
  std::array<int, 3> pbc{{true, true, true}};
  cell << 6.19, 2.41, 0.21, 0.00, 6.15, 1.02, 0.00, 0.00, 7.31;

  // clang-format off
  positions << 3.689540159937393, 5.123016813620886, 1.994119731169116,
      6.818437242389163, 2.630056617829216, 6.182500355729062,
      2.114977334498767, 6.697579639059512, 1.392155450018263,
      7.420401523540017, 2.432242071439904, 6.380314902118375,
      1.112656394115962, 7.699900579442317, 3.569715877854675,
      5.242841095703604, 3.122826344932127, 5.689730628626151,
      3.248684682453303, 5.563872291104976, 2.608353462112637,
      6.204203511445642, 5.035681855581504, 2.134827911489532,
      0.946910011088814, 6.223599755982222, 4.168634519120968,
      3.001875247950068, 1.980327734683430, 5.190182032387606,
      2.943861424421339, 4.226648342649697, 5.457161501166098,
      1.713348265904937, 1.501663178733906, 5.668846588337130,
      5.208365510425203, 1.962144256645833, 2.728127406527150,
      4.442382360543885, 2.839975217222644, 4.330534549848392,
      0.744216089807768, 6.426293677263268, 4.643695520786083,
      2.662204050783991, 1.250682335857938, 6.055217235712136,
      0.860905287815103, 6.444994283754972, 4.536108843695142,
      2.769790727874932, 5.609177455068640, 1.696722116501434,
      6.703053268421970, 0.602846303148105, 3.487609972580834,
      3.818289598989240, 1.436734374347541, 5.869165197222533,
      1.054504320562138, 6.251395251007936, 3.998423858825871,
      3.307475712744203, 5.323662899811682, 1.982236671758393;
  // clang-format on
  positions.transposeInPlace();
  atom_types << 20, 20, 24, 24, 15, 15, 15, 15, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8;

  using PBC_t = Eigen::Map<Eigen::Matrix<int, 3, 1>>;
  // manager->update(positions, atom_types, cell, PBC_t{pbc.data()});
  auto manager{rascal::make_structure_manager<Manager_t>()};
  double cutoff{2.5};
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

  // std::string filename{"crystal_structure.json"};
  // std::string filename{"alanine-X.json"};
  std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  //"reference_data/CaCrP2O7_mvc-11955_symmetrized_.json",

  std::ifstream reader(filename.c_str());
  std::cout << reader.is_open() << std::endl;
  std::string str((std::istreambuf_iterator<char>(reader)),
                 std::istreambuf_iterator<char>());
  std::cout << str << std::endl;
  std::cout << "######################" << std::endl;

  json j;
  reader.seekg(0);
  reader >> j;

  std::cout << j.dump(2) << std::endl;
  std::cout << "Reading structure " << filename << std::endl;

  // `manager` is the data object which reads, stores and gives access to atom
  // positions, types. It also provides iteration over all atom.
//  manager.update(filename);

  // `pair_manager` is constructed with the `manager` and a `cutoff`.
  auto pair_manager{rascal::make_adapted_manager<rascal::AdaptorNeighbourList>(manager, cutoff, true)};
  // By invoking the `.update()` method, a neighbour list is built.
//  pair_manager->update();

  // `strict_manager` is constructed with a `pair_manager`.
  auto strict_manager{
        rascal::make_adapted_manager<rascal::AdaptorStrict>(pair_manager, cutoff)};
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
    std::cout << "atom " << atom.get_atom_index() << " global index "
              << atom.get_global_index() << std::endl;
  }

  // `pair_manager` provides iteration over atoms and pairs
  for (auto atom : pair_manager) {
    for (auto pair : atom) {
      std::cout << "pair (" << atom.get_atom_index() << ", "
                << pair.get_atom_index() << " ) global index "
                << pair.get_global_index() << std::endl;
    }
  }

  // `strict_manager` provides iteration over atoms and strict pairs
  for (auto atom : strict_manager) {
    for (auto pair : atom) {
      std::cout << "strict pair (" << atom.get_atom_index() << ", "
                << pair.get_atom_index() << ") global index "
                << pair.get_global_index() << std::endl;
    }
  }

  // `triplet_manager` provides iteration over atoms, strict pairs and strict
  // triplets
  for (auto atom : triplet_manager) {
    for (auto pair : atom) {
      for (auto triplet : pair) {
        std::cout << "triplet (" << atom.get_atom_index() << ", "
                  << pair.get_atom_index() << ", " << triplet.get_atom_index()
                  << ") global index " << triplet.get_global_index()
                  << std::endl;
      }
    }
  }
}
