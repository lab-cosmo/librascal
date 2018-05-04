/**
 * file   test_properties.cc
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  tests for cluster-related properties
 *
 * @section LICENSE
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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
#include "test_neighbourhood.hh"
#include "neighbourhood_managers/property.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE (Property_tests);
  struct PropertyFixture: public ManagerFixture {
    using Manager_t = typename ManagerFixture::Manager_t;
    using PairScalarProperty_t = Property<Manager_t, double, 2>;
    using AtomVectorProperty_t = Property<Manager_t, double, 1, 3>;


    PropertyFixture()
      :ManagerFixture{}, pair_property{this->manager},
       atom_property{this->manager}
    {}

    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
  };

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test, PropertyFixture) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(fill_test, PropertyFixture) {
    pair_property.resize();
    atom_property.resize();
    int pair_property_counter{};
    for (auto atom: manager) {
      atom_property[atom] = atom.get_x();
      for (auto pair: atom) {
        pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom: manager) {
      auto error = (atom_property[atom] - atom.get_x()).norm();
      BOOST_CHECK_EQUAL(error, 0);
      for (auto pair: atom) {
        BOOST_CHECK_EQUAL(pair_property[pair], ++pair_property_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(compute_distances, PropertyFixture) {
    pair_property.resize();

    for (auto atom: manager) {
      for (auto pair: atom) {
        pair_property[pair] = (atom.get_x() - pair.get_x()).norm();
      }
    }

    for (auto atom: manager) {
      for (auto pair: atom) {
        auto error = pair_property[pair] - 1;
        BOOST_CHECK_LE(error, 1e-12);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END ();

}  // rascal
