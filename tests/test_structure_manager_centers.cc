/**
 * @file     test_structure_manager_centers.cc
 *
 * @author  Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date    09 Aug 2018
 *
 * @brief   Test the implementation of the Order=1 center manager
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

namespace rascal {

  BOOST_AUTO_TEST_SUITE(ManagerCentersTests);
  /* ---------------------------------------------------------------------- */
  // checking the constructor
  BOOST_FIXTURE_TEST_CASE(manager_centers_constructor_test,
                          ManagerFixture<StructureManagerCenters>) {}

  /* ---------------------------------------------------------------------- */
  /**
   * checking the iteration with and without the center mask
   *
   * test that the positions are properly iterated over
   */
  BOOST_FIXTURE_TEST_CASE(iterator_test,
                          ManagerFixture<StructureManagerCenters>) {
    constexpr bool verbose{false};
    int i_manager{0};
    for (auto & manager : this->managers) {
      auto & structure = this->structures[i_manager];
      auto & center_atoms_mask = structure.center_atoms_mask;
      BOOST_CHECK_EQUAL(manager->get_size(), center_atoms_mask.count());
      BOOST_CHECK_EQUAL(manager->get_nb_clusters(1), center_atoms_mask.count());

      int atom_counter{0};

      auto & positions = structure.positions;
      for (auto atom_cluster : manager) {
        BOOST_CHECK_EQUAL(atom_counter, atom_cluster.get_index());
        auto index = manager->get_atom_index(atom_cluster);
        auto error =
            (positions.col(index) - atom_cluster.get_position()).norm();
        BOOST_CHECK_LE(error, 1e-14);

        if (verbose) {
          std::cout << "atom (" << atom_cluster.back()
                    << "), atom counter = " << atom_counter << std::endl;
        }
        ++atom_counter;
      }
      ++i_manager;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Check that the atom in unit cell check is properly done
   */
  BOOST_FIXTURE_TEST_CASE(atom_bound_test,
                          ManagerFixture<StructureManagerCenters>) {
    int i_manager{0};
    for (auto & manager : this->managers) {
      auto & structure = this->structures[i_manager];
      structure.positions.col(0) << -0.5, -0.5, -0.5;
      BOOST_CHECK_THROW(manager->update(structure), std::runtime_error);
      ++i_manager;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test for checking StructureManagerCenters specific interface, ``manager``
   * is the name of the manager object. it checks specifically, if the data
   * which is read from a json file (or in the case
   */
  BOOST_FIXTURE_TEST_CASE(simple_cubic_9_neighbour_list,
                          ManagerFixtureFile<StructureManagerCenters>) {
    constexpr bool verbose{false};
    for (auto & manager : this->managers) {
      if (verbose)
        std::cout << "StructureManagerCenters interface" << std::endl;

      auto && dim{manager->dim()};
      if (verbose) {
        std::cout << "dimension: " << dim << std::endl;
      }

      auto && cell{manager->get_cell()};
      if (verbose) {
        std::cout << "cell:\n" << cell << std::endl;
      }

      auto && atom_types{manager->get_atom_types()};
      if (verbose) {
        std::cout << "atom types:\n" << atom_types << std::endl;
      }

      auto && positions{manager->get_positions()};
      if (verbose) {
        std::cout << "atom positions:\n" << positions << std::endl;
      }

      auto && periodicity{manager->get_periodic_boundary_conditions()};
      if (verbose) {
        std::cout << "periodicity (x,y,z):\n" << periodicity << std::endl;
      }

      auto center_atoms_mask{manager->get_center_atoms_mask()};
      if (verbose) {
        std::cout << "which atom is a centeral atom:\n"
                  << center_atoms_mask << std::endl;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  // checking update
  BOOST_FIXTURE_TEST_CASE(manager_update_test,
                          ManagerFixture<StructureManagerCenters>) {
    int i_manager{0};
    for (auto & manager : this->managers) {
      auto natoms = manager->size();
      auto natoms2 = manager->get_size();
      auto natoms3 = manager->get_nb_clusters(1);
      BOOST_CHECK_EQUAL(natoms, natoms2);
      BOOST_CHECK_EQUAL(natoms, natoms3);
      auto & structure = this->structures[i_manager];
      auto & positions = structure.positions;
      auto & atom_types = structure.atom_types;

      for (auto atom : manager) {
        auto index = manager->get_atom_index(atom);
        auto type = atom.get_atom_type();
        BOOST_CHECK_EQUAL(type, atom_types[index]);

        auto cluster_size = manager->get_cluster_size(atom);
        BOOST_CHECK_EQUAL(cluster_size, 1);

        auto pos = atom.get_position();
        auto pos_reference = positions.col(index);
        auto position_error = (pos - pos_reference).norm();
        BOOST_CHECK(position_error < tol / 100);
      }
      ++i_manager;
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
