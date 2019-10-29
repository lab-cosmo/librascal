/**
 * @file   test_key_standardisation.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 Apr 2019
 *
 * @brief  test configuration for key standardisation
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

#include "utils/key_standardisation.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {
  /**
   * Overload of stdout stream for tuple standardisation for debugging purposes
   */
  template <typename T, size_t Order>
  std::ostream & operator<<(std::ostream & os,
                            const KeyStandardisation<T, Order> & index) {
    os << "(";
    for (size_t i = 0; i < Order - 1; ++i) {
      os << index[i] << ", ";
    }
    os << index.back() << ")";
    return os;
  }

  BOOST_AUTO_TEST_SUITE(key_standardisation_test);
  /* ---------------------------------------------------------------------- */
  /**
   * Fixture to provide the an integer array used for the construction of a
   * `KeyStandardisation` type based on a given `Order`.
   */
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

  /**
   * Fixture to provide template `Order`/`MaxOrder` combinations for a
   * `KeyStandardisation` type.
   */
  template <size_t Order, size_t MaxOrder>
  struct KeyFixture {
    KeyFixture() : standard_key{indices_fixture.indices} {}
    ~KeyFixture() = default;

    IndexFixture<Order> indices_fixture{};
    KeyStandardisation<int, MaxOrder> standard_key;
    static size_t get_order() { return Order; }
  };

  /**
   * List for testing different combinations of Order and MaxOrder to check the
   * reinterpret_cast with the bracket operator.
   */
  using KeyFixtures =
      boost::mpl::list<KeyFixture<1, 2>, KeyFixture<2, 2>, KeyFixture<1, 3>,
                       KeyFixture<2, 3>, KeyFixture<3, 3>, KeyFixture<7, 8>>;

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, KeyFixtures, Fix) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(iterator_test, Fix, KeyFixtures, Fix) {
    for (size_t i{0}; i < Fix::get_order(); ++i) {
      auto ref_val{Fix::indices_fixture.indices[i]};
      auto key_val{Fix::standard_key[i]};
      BOOST_CHECK_EQUAL(ref_val, key_val);
    }
  }

  BOOST_AUTO_TEST_SUITE_END()
}  // namespace rascal
