/**
 * file   test_structure_manager.cc
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
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_structure.hh"

#include <iostream>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerJsonFileTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_constructor_test,
                          ManagerFixtureFile<StructureManagerCenters>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test,
                          ManagerFixtureFile<StructureManagerCenters>) {
    // Reference values
    constexpr int natoms{9};

    BOOST_CHECK_EQUAL(manager.get_size(), natoms);
    BOOST_CHECK_EQUAL(manager.get_nb_clusters(1), natoms);

    int atom_counter{};
    //    constexpr bool verbose{false};

    for (auto atom_cluster : manager) {
      BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_global_index());
      ++atom_counter;
    }
    BOOST_CHECK_EQUAL(atom_counter, natoms);
  }
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(test_get_cell,
                          ManagerFixtureFile<StructureManagerCenters>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
