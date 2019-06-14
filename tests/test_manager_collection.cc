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

  using fixtures_test = boost::mpl::list<ManagerNLCollectionFixture>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(load_multiple_structures_test, Fix,
                                   fixtures_test, Fix) {
    auto& collection = Fix::collection;
    auto& filename = Fix::filename;

    collection.add_structures(filename);

    std::cout << collection.size() << std::endl;

  }
  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal

