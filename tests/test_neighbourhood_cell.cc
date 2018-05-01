/**
 * file   lattice.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of the neighbourhood manager for lammps
 *        neighbourhood lists
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
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
#include "test_neighbourhood_cell.hh"


namespace proteus {

  BOOST_AUTO_TEST_SUITE(ManagerCellTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_cell_constructor_test, ManagerFixture_cell){}
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_test, ManagerFixture_cell) {
    
    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{false};
    for (auto atom_cluster: manager) {
      BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_index());
      ++atom_counter;
      /*
      for (auto pair_cluster: atom_cluster) {
        auto pair_offset{pair_cluster.get_global_index()};
        if (verbose) {
          std::cout << "pair (" << atom_cluster.get_atoms().back().get_index()
                    << ", " << pair_cluster.get_atoms().back().get_index()
                    << "), pair_counter = " << pair_counter
                    << ", pair_offset = " << pair_offset << std::endl;
        }

        BOOST_CHECK_EQUAL(pair_counter, pair_offset);
        ++pair_counter;

      }*/
    }
    
  }
  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();
}  // proteus
