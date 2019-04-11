/**
 * file   test_tuple_standardisation.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 Apr 2019
 *
 * @brief  test configuration for tuple standardisation
 *
 * Copyright © 2019 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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
 * along with this software; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "utils/tuple_standardisation.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(tuple_standardisation_test);

  template <size_t Order>
  struct IndexFixture {
    IndexFixture() {
      for (size_t i{0}; i < Order; ++i) {
        indices[i] = i;
      }
    }
    ~IndexFixture() = default;
    std::array<int, Order> indices{};
  };

  template <size_t Order>
  struct TupleFixture {
    TupleFixture() : standard_tuple{indices_fixture.indices} {}
    ~TupleFixture() = default;

    IndexFixture<Order> indices_fixture{};
    TupleStandardisation<int, Order> standard_tuple;
  };

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test, TupleFixture<5>) {}

  /* ---------------------------------------------------------------------- */

  BOOST_FIXTURE_TEST_CASE(iterator_test, TupleFixture<5>) {
    std::cout << "indices " << std::endl;
    for (size_t i{0}; i < 5; ++i) {
      std::cout << indices_fixture.indices[i] << " ";
    }
    // for (auto && i : standard_tuple) {
    //   std::cout << i << std::endl;
    // }

    std::cout << std::endl;
    std::cout << "StandardTuple " << standard_tuple << std::endl;
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // namespace rascal
