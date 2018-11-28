/**
 * file   test_structure_manager_chain.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 Aug 2018
 *
 * @brief  tests for the class NeighbourhoodManagerJson
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_structure.hh"

#include <iostream>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerJsonTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_json_constructor_test,
                          ManagerFixture<StructureManagerJson>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test, ManagerFixture<StructureManagerJson>) {
    // Reference values
    constexpr int natoms{9};

    BOOST_CHECK_EQUAL(manager_json.get_size(), natoms);
    BOOST_CHECK_EQUAL(manager_json.get_nb_clusters(1), natoms);

    int atom_counter{};
    //    constexpr bool verbose{false};

    for (auto atom_cluster : manager_json) {
      BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_global_index());
      ++atom_counter;
    }
    BOOST_CHECK_EQUAL(atom_counter, natoms);
  }
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(test_get_cell, ManagerFixture<StructureManagerJson>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
