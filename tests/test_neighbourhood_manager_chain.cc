/**
 * file   test_neighbourhood_manager_chain.cc
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
#include "test_neighbourhood.hh"

#include <iostream>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerChainTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_chain_constructor_test,
			  ManagerFixture_chain) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test, ManagerFixture_chain) {
    // Reference values
    constexpr int natoms{9};
    constexpr int npairs{36};

    // std::cout << "Chain: test " << manager_chain.get_nb_clusters(2) << "\n";
    // Check number of atoms
    BOOST_CHECK_EQUAL(manager_chain.get_size(), natoms);
    // Check cluster Level 1, atoms
    BOOST_CHECK_EQUAL(manager_chain.get_nb_clusters(1), natoms);
    // Check cluster Level 2, pairs
    BOOST_CHECK_EQUAL(manager_chain.get_nb_clusters(2), npairs);

    int atom_counter{0};
    int pair_counter{0};
    constexpr bool verbose{true};

    for (auto atom_cluster : manager_chain) {
      BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_global_index());
      ++atom_counter;
      for (auto pair_cluster : atom_cluster) {
      	auto pair_offset{pair_cluster.get_global_index()};
      	++pair_counter;
      	if (verbose) {
      	  std::cout
	    << "pair (" << atom_cluster.get_atoms().back().get_index()
	    << ", " << pair_cluster.get_atoms().back().get_index()
	    << "), pair_counter = " << pair_counter
	    << ", pair_offset = " << pair_offset << std::endl;
        }
      }
    }
    BOOST_CHECK_EQUAL(atom_counter, natoms);
    BOOST_CHECK_EQUAL(pair_counter, npairs);
  }
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test, ManagerFixture_chain) {

  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
