/**
 * file   test_structure_manager_lammps.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  tests for the class `NeighbourhoodManagerLammps
 *
 * @section LICENSE
 *
 * Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#include "tests.hh"
#include "test_structure.hh"

#include <iostream>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerLammps>) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test,
                          ManagerFixture<StructureManagerLammps>) {
    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{false};

    for (auto atom : manager) {
      BOOST_CHECK_EQUAL(atom.get_global_index(), atom_counter);
      BOOST_CHECK_EQUAL(atom.get_atom_type(),
                        manager->get_atom_type(atom.get_atom_tag()));
      ++atom_counter;

      for (auto pair : atom) {
        auto pair_offset{pair.get_global_index()};
        if (verbose) {
          std::cout << "pair (" << atom.get_atom_tag() << ", "
                    << pair.get_atom_tag()
                    << "), pair_counter = " << pair_counter
                    << ", pair_offset = " << pair_offset << std::endl;
        }

        BOOST_CHECK_EQUAL(pair_offset, pair_counter);
        ++pair_counter;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(interface_test,
                          ManagerFixture<StructureManagerLammps>) {
    auto natoms = manager->size();
    auto natoms2 = manager->get_size();
    BOOST_CHECK_EQUAL(natoms, natoms2);

    for (auto atom : manager) {
      auto atom_tag = atom.get_atom_tag();
      auto atom_type = atom.get_atom_type();
      BOOST_CHECK_EQUAL(atom_type, manager->get_atom_type(atom_tag));
      // auto index_size = manager->get_cluster_size(index);
      // auto cluster_size = manager->get_cluster_size(atom);
      // BOOST_CHECK_EQUAL(index_size, cluster_size);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
