/**
 * file test_manager_collection.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 June 2019
 *
 * @brief test the structure manager collection
 *
 * @section LICENSE
 *
 * Copyright  2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_manager_collection.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(manager_collection_test);

  using fixtures_test =
      boost::mpl::list<CollectionFixture<StrictNLCollectionFixture>>;

  /**
   * Test loading structures from file in ubjson binary format
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(load_multiple_structures_test, Fix,
                                   fixtures_test, Fix) {
    auto & collections = Fix::collections;
    auto & filename = Fix::filename;
    auto & start = Fix::start;
    auto & length = Fix::length;

    for (auto & collection : collections) {
      collection.add_structures(filename, start, length);
      BOOST_CHECK_EQUAL(collection.size(), length);
    }
  }

  /**
   * Test the iteration over the manager and its functionalities
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(iterate_structures_test, Fix, fixtures_test,
                                   Fix) {
    auto & collections = Fix::collections;
    auto & filename = Fix::filename;
    auto & start = Fix::start;
    auto & length = Fix::length;

    for (auto & collection : collections) {
      // add some integer constant to not start at position zero
      collection.add_structures(filename, start + 30, length + 5);

      for (auto & manager : collection) {
        for (auto atom : manager) {
          BOOST_CHECK_EQUAL(atom.get_position().size(), 3);
        }
      }
    }
  }

  /**
   * Test adding structures with json hyper parameter format
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(add_structures_test, Fix, fixtures_test,
                                   Fix) {
    auto & collections = Fix::collections;
    auto & structures = Fix::structures;
    for (auto & collection : collections) {
      collection.add_structures(structures);
      BOOST_CHECK_EQUAL(collection.size(), structures.size());
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
