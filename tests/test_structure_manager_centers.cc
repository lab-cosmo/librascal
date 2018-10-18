/**
 * file     test_structure_manager_centers.cc
 *
 * @author  Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date    09 Aug 2018
 *
 * @brief   Test the implementation of the Order=1 center manager
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerCentersTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(manager_centers_constructor_test,
                          ManagerFixture<StructureManagerCenters>){}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test,
                          ManagerFixture<StructureManagerCenters>) {

    BOOST_CHECK_EQUAL(manager.get_size(), numbers.size());
    BOOST_CHECK_EQUAL(manager.get_nb_clusters(1), numbers.size());

    int atom_counter{};
    constexpr bool verbose{false};

    for (auto atom_cluster: manager) {
      BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_index());
      ++atom_counter;
      for (int ii{3};ii<3;++ii){
        BOOST_CHECK_EQUAL(positions(ii,atom_cluster.get_index()),
                          atom_cluster.get_position()[ii]);
      }
      if (verbose) std::cout << "atom (" << atom_cluster.back()
                             << "), atom counter = "
                             << atom_counter
                             << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(interface_test,
                          ManagerFixture<StructureManagerCenters>) {

    auto natoms = manager.size();
    auto natoms2 = manager.get_size();
    auto natoms3 = manager.get_nb_clusters(1);
    BOOST_CHECK_EQUAL(natoms, natoms2);
    BOOST_CHECK_EQUAL(natoms, natoms3);

    for (auto atom : manager) {
      auto index = atom.get_atom_index();
      auto type = atom.get_atom_type();
      BOOST_CHECK_EQUAL(type, numbers[index]);

      auto cluster_size = manager.get_cluster_size(atom);
      BOOST_CHECK_EQUAL(cluster_size, 1);

      auto pos = atom.get_position();
      auto pos_reference = positions.col(index);
      auto position_error = (pos - pos_reference).norm();
      BOOST_CHECK(position_error < tol / 100);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();
}  // rascal
