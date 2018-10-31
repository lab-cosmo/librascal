/**
 * file     test_structure_managers_interface.cc
 *
 * @author  Michele Ceriotti <michele.ceriotti@epfl.ch>
 *
 * @date    22 Oct 2018
 *
 * @brief   A template test for the basic interface of structure managers
 *
 * Copyright © 2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
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

  BOOST_AUTO_TEST_SUITE(structure_managers_interface);

  // gets a list of fixtures for all the different possible structure managers
  using fixtures = boost::mpl::list<
    ManagerFixture<StructureManagerCenters>,
    ManagerFixture<StructureManagerLammps>,
    ManagerFixtureFile<StructureManagerCenters>,
    // not templated single manager fixtures
    ManagerFixtureSimple>;

  /* ---------------------------------------------------------------------- */
  // just checks that the structure managers can be constructed
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_constructor_test, Fix, fixtures,
                                   Fix) { }

  /* ---------------------------------------------------------------------- */
  // the size of the manager should correspond to the number of size 1 clusters
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_size_consistency, Fix, fixtures,
                                   Fix) {
    auto & manager = Fix::manager;
    BOOST_CHECK_EQUAL(manager.size(), manager.nb_clusters(1));
  }

  /* ---------------------------------------------------------------------- */
  // loops over the centers in the manager making sure positions are consistent
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_atom_indexing, Fix, fixtures,
                                   Fix) {
    auto & manager = Fix::manager;
    for (auto atom : manager) {
      // checks get_atom_index exists
      auto index = atom.get_atom_index();

      // checks get_atom_type exists
      auto type = atom.get_atom_type();

      // cluster size should be 1!
      BOOST_CHECK_EQUAL(type, manager.atom_type(index));

      // checks that multiple ways of accessing positions are equivalent
      auto position_error = (atom.get_position() -
                             manager.position(index)).norm();

      BOOST_CHECK(position_error < tol / 100);

      position_error = (atom.get_position() -
                        manager.position(atom.back())).norm();
      BOOST_CHECK(position_error < tol / 100);
    }
  }

  /* ---------------------------------------------------------------------- */
  // loops over the centers in the manager making global atom indices are
  // contiguous
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_atom_global_indexing, Fix,
                                   fixtures, Fix) {
    auto & manager = Fix::manager;

    auto index_reference{0};

    for (auto atom : manager) {
      // checks get_atom_index exists
      auto index = atom.get_global_index();
      BOOST_CHECK_EQUAL(index_reference, index);
      index_reference++;
    }
  }


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();
}  // rascal
