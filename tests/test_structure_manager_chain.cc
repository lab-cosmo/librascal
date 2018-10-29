/**
 * file   test_structure_manager.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   06 Jun 2018
 *
 * @brief  tests for the class NeighbourhoodManagerChain
 *
 * @section LICENSE
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_structure.hh"

#include <iostream>

namespace rascal {

  // BOOST_AUTO_TEST_SUITE(ManagerChainTests);
  // /* ---------------------------------------------------------------------- */
  // BOOST_FIXTURE_TEST_CASE(manager_constructor_test,
  //                         ManagerFixture<StructureManagerChain>) {
  // }

  // /* ---------------------------------------------------------------------- */
  // BOOST_FIXTURE_TEST_CASE(iterator_test,
  //                         ManagerFixture<StructureManagerChain>) {
  //   // Reference values
  //   constexpr static int natoms{9};
  //   constexpr static int npairs{36};

  //   BOOST_CHECK_EQUAL(manager.get_size(), natoms);
  //   BOOST_CHECK_EQUAL(manager.get_nb_clusters(1), natoms);
  //   BOOST_CHECK_EQUAL(manager.get_nb_clusters(2), npairs);

  //   int atom_counter{};
  //   int pair_counter{};
  //   constexpr bool verbose{false};

  //   for (auto atom_cluster : manager) {
  //     BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_global_index());
  //     ++atom_counter;
  //     for (auto pair_cluster : atom_cluster) {
  //       auto pair_offset{pair_cluster.get_global_index()};
  //       if (verbose) {
  //         std::cout
  //           << "pair (" << atom_cluster.back()
  //           << ", " << pair_cluster.back()
  //           << "), pair_counter = " << pair_counter
  //           << ", pair_offset = " << pair_offset << std::endl;
  //       }
  //       BOOST_CHECK_EQUAL(pair_counter, pair_offset);
  //       ++pair_counter;
  //     }
  //   }
  //   BOOST_CHECK_EQUAL(atom_counter, natoms);
  //   BOOST_CHECK_EQUAL(pair_counter, npairs);
  // }
  // /* ---------------------------------------------------------------------- */
  // // BOOST_FIXTURE_TEST_CASE(neighbourlist_test, ManagerFixture<StructureManagerChain>) {
  // //   std::vector<std::vector<int>> neighbours(9);
  // //   neighbours[0] = {1, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[1] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[2] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[3] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[4] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[5] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[6] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[7] = {0, 5, 3, 7, 2, 6, 4, 8};
  // //   neighbours[8] = {0, 5, 3, 7, 2, 6, 4, 8};

  // //   for (size_t i=0; i < manager.get_size(); ++i) {
  // //     neighs = manager.
  // //   }


  // // }
  // /* ---------------------------------------------------------------------- */
  // BOOST_FIXTURE_TEST_CASE(test_get_cell, ManagerFixture<StructureManagerChain>) {
  // }

  // /* ---------------------------------------------------------------------- */
  // BOOST_AUTO_TEST_SUITE_END();

}  // rascal
