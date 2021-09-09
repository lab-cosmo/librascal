/**
 * @file   test_permutation.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   09 Sep 2021
 *
 * @brief  tests for the compile-time permutation handlers
 *
 * Copyright © 2021 Till Junge
 *
 * LibRascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * LibRascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LibRascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it
 * with proprietary FFT implementations or numerical libraries, containing parts
 * covered by the terms of those libraries' licenses, the licensors of this
 * Program grant you additional permission to convey the resulting work.
 *
 */

#include <array>

#include "rascal/utils/permutation.hh"

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(permutation_tests);

  template <size_t I, size_t J, size_t K, bool a, bool b, bool c>
  struct PermutationTriplets {
    using type = Permutation<3, I, J, K>;
    const std::array<bool,3> pair_inversions { a, b, c };
  };

  template <size_t I, size_t J, bool inv>
  struct PermutationPairs {
    using type = Permutation<2, I, J>;
    const std::array<bool, 1> pair_inversions { inv };
  };

  using Permutations =
      boost::mpl::list<PermutationTriplets<0, 1, 2, false, false, true>,
                       PermutationTriplets<1, 2, 0, false, true, false>,
                       PermutationTriplets<2, 0, 1, true, false, false>,
                       PermutationTriplets<0, 2, 1, false, true, true>,
                       PermutationTriplets<2, 1, 0, true, true, false>,
                       PermutationTriplets<1, 0, 2, true, false, true>,
                       PermutationPairs<0, 1, false>,
                       PermutationPairs<1, 0, true>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(test_pair_inversion, Fix, Permutations, Fix) {
    const auto inversion { Fix::type::pair_inversion() };

    for (size_t i{0}; i < inversion.size(); ++i) {
      BOOST_CHECK_EQUAL(inversion[i], Fix::pair_inversions[i]);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
