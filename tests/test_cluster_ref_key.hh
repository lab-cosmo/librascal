/**
 * @file   test_cluster_ref_key.hh
 *
 * @author Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   10 Oct 2019
 *
 * @brief  tests for the cluster ref key layer computation utilities
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef TESTS_TEST_CLUSTER_REF_KEY_HH_
#define TESTS_TEST_CLUSTER_REF_KEY_HH_

#include "rascal/structure_managers/cluster_ref_key.hh"

#include <tuple>
#include <vector>

namespace rascal {
  struct LayerFixture {
    // {1}
    using order_one_seq = std::index_sequence<1>;
    // {1,1}
    using order_two_seq = std::index_sequence<1, 1>;
    // {1,1,0}
    using order_three_seq = std::index_sequence<1, 1, 0>;
  };

  /**
   * Helper function which prints the index sequence `Ints` defined in the
   * std::index_sequence type of the parameter `sequence` by recursively
   * printing the head element.
   *
   * @param sequence sequence to be printed
   */
  template <size_t... Ints>
  void print_index_sequence(std::index_sequence<Ints...> /*sequence*/) {
    std::vector<size_t> vec = {Ints...};
    std::cout << "( ";
    for (size_t ele : vec) {
      std::cout << ele << " ";
    }
    std::cout << ")" << std::endl;
  }

  /**
   * Helper function which prints the index sequence defined in the
   * std::index_sequence type of the parameter `sequence` up to template
   * parameter `Length`
   *
   * @tparam Length
   * @param sequence sequence to be printed
   */
  template <size_t Length, size_t... Ints>
  static void print_index_sequence(std::index_sequence<Ints...> /*sequence*/) {
    std::cout << std::array<size_t, Length>{Ints...} << std::endl;
  }
}  // namespace rascal
#endif  // TESTS_TEST_CLUSTER_REF_KEY_HH_
