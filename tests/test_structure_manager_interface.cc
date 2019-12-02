/**
 * @file     test_structure_managers_interface.cc
 *
 * @author  Michele Ceriotti <michele.ceriotti@epfl.ch>
 *
 * @date    22 Oct 2018
 *
 * @brief   A template test for the basic interface of structure managers
 *
 * Copyright  2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

constexpr double TOLERANCE = 1e-14;

namespace rascal {

  BOOST_AUTO_TEST_SUITE(structure_managers_interface);

  /* ---------------------------------------------------------------------- */
  // a list of fixtures for all the different possible structure managers
  // and test the atom (Order=1) related interface.
  using fixtures = boost::mpl::list<ManagerFixture<StructureManagerCenters>,
                                    ManagerFixture<StructureManagerLammps>,
                                    ManagerFixtureFile<StructureManagerCenters>,
                                    // not templated single manager fixtures
                                    ManagerFixtureSimple>;

  /* ---------------------------------------------------------------------- */
  // a list of fixtures for pair managers for testing pair related interface
  using pair_fixtures =
      boost::mpl::list<ManagerFixture<StructureManagerLammps>>;

  /* ---------------------------------------------------------------------- */
  // just checks that the structure managers can be constructed
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_constructor_test, Fix, fixtures,
                                   Fix) {}

  /* ---------------------------------------------------------------------- */
  // the size of the manager should correspond to the number of size 1 clusters
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_size_consistency, Fix, fixtures,
                                   Fix) {
    auto & manager = Fix::manager;
    BOOST_CHECK_EQUAL(manager->size(), manager->nb_clusters(1));
  }

  /* ---------------------------------------------------------------------- */
  // loops over the centers in the manager making sure positions are consistent
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_atom_taging, Fix, fixtures, Fix) {
    auto & manager = Fix::manager;
    for (auto atom : manager) {
      // checks get_atom_tag exists
      auto index = atom.get_atom_tag();

      // checks get_atom_type exists
      auto type = atom.get_atom_type();

      // cluster size should be 1!
      BOOST_CHECK_EQUAL(type, manager->atom_type(index));

      // checks that multiple ways of accessing positions are equivalent
      auto position_error =
          (atom.get_position() - manager->position(index)).norm();

      BOOST_CHECK(position_error < TOLERANCE);

      position_error =
          (atom.get_position() - manager->position(atom.back())).norm();
      BOOST_CHECK(position_error < TOLERANCE);
    }
  }

  /* ---------------------------------------------------------------------- */
  // loops over the centers in the manager making global atom tags are
  // contiguous
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_atom_global_indexing, Fix,
                                   fixtures, Fix) {
    auto & manager = Fix::manager;

    auto index_reference{0};

    for (auto atom : manager) {
      // checks get_atom_tag exists
      auto index = atom.get_global_index();
      BOOST_CHECK_EQUAL(index_reference, index);
      index_reference++;
    }
  }

  /* ---------------------------------------------------------------------- */
  // test template for pairs
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_pair_global_indexing, Fix,
                                   pair_fixtures, Fix) {
    auto & manager = Fix::manager;
    auto pair_reference{0};

    for (auto atom : manager) {
      for (auto pair : atom) {
        auto global_index = pair.get_global_index();
        BOOST_CHECK_EQUAL(pair_reference, global_index);
        pair_reference++;

        // check index access
        auto index = pair.get_atom_tag();

        // check atom type access
        auto type = pair.get_atom_type();

        BOOST_CHECK_EQUAL(type, manager->atom_type(index));

        // check positions
        auto position_error =
            (pair.get_position() - manager->position(index)).norm();
        BOOST_CHECK(position_error < TOLERANCE);

        position_error =
            (pair.get_position() - manager->position(pair.back())).norm();
        BOOST_CHECK(position_error < TOLERANCE);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
