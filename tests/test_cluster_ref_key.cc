/**
 * @file   test_cluster_ref_key.cc
 *
 * @author Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   10 Oct 2019
 *
 * @brief  tests for the cluster ref key layer computation utilities
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_cluster_ref_key.hh"

namespace rascal {
  BOOST_AUTO_TEST_SUITE(test_cluster_ref_key);

  /**
   * tests the get_min_layer function
   */
  BOOST_AUTO_TEST_CASE(test_get_min_layer) {
    // LayersByOrder using only one root manager
    const size_t active_layer_root_manager{
        get_min_layer<1>(typename std::index_sequence<0>{})};
    BOOST_CHECK_EQUAL(0, active_layer_root_manager);
    // LayersByOrder for `AdaptorFullNeighbourList`
    const size_t active_layer_order_one_adaptor{
        get_min_layer<1>(typename std::index_sequence<1, 0>{})};
    BOOST_CHECK_EQUAL(1, active_layer_order_one_adaptor);
    // LayersByOrder for stacking an `AdaptorMaxOrder` on an `AdaptorStrict`
    const size_t active_layer_increase_max_order{
        get_min_layer<3>(typename std::index_sequence<1, 1, 0>{})};
    BOOST_CHECK_EQUAL(0, active_layer_increase_max_order);
  }

  /**
   * tests the get_layer function
   */
  BOOST_FIXTURE_TEST_CASE(test_get_layer, LayerFixture) {
    BOOST_CHECK(get_layer<1>(order_three_seq{}) == 1);
    BOOST_CHECK(get_layer<2>(order_three_seq{}) == 1);
    BOOST_CHECK(get_layer<3>(order_three_seq{}) == 0);
  }

  /**
   * tests the index_sequence_to_array function
   */
  BOOST_FIXTURE_TEST_CASE(test_index_sequence_to_array, LayerFixture) {
    bool verbose{false};

    std::array<size_t, 1> arr_order_one =
        index_sequence_to_array(order_one_seq{});
    std::array<size_t, 1> arr_order_one_ref = {1};
    BOOST_CHECK(arr_order_one == arr_order_one_ref);

    std::array<size_t, 2> arr_order_two =
        index_sequence_to_array(order_two_seq{});
    std::array<size_t, 2> arr_order_two_ref = {1, 1};
    BOOST_CHECK(arr_order_two == arr_order_two_ref);

    std::array<size_t, 3> arr_order_three =
        index_sequence_to_array(order_three_seq{});
    std::array<size_t, 3> arr_order_three_ref = {1, 1, 0};
    BOOST_CHECK(arr_order_three == arr_order_three_ref);
    if (verbose) {
      print_index_sequence(std::make_index_sequence<5>{});
    }
  }

  /**
   * tests the LayerExtender struct
   */
  BOOST_AUTO_TEST_CASE(test_layer_extender) {
    bool verbose{false};

    using LayersUpToOrderTwo =
        typename LayerExtender<std::index_sequence<1>>::type;
    using LayersUpToOrderThree =
        typename LayerExtender<LayersUpToOrderTwo>::type;

    if (verbose) {
      print_index_sequence(LayersUpToOrderTwo{});
    }
    BOOST_CHECK(typeid(std::index_sequence<1, 0>) ==
                typeid(LayersUpToOrderTwo));

    if (verbose) {
      print_index_sequence(LayersUpToOrderThree{});
    }
    BOOST_CHECK(typeid(std::index_sequence<1, 0, 0>) ==
                typeid(LayersUpToOrderThree));
  }

  /**
   * tests the LayerIncreaser struct
   */
  BOOST_AUTO_TEST_CASE(test_layer_increaser) {
    bool verbose{false};
    using LayersUpToOrderThree =
        typename LayerIncreaser<std::index_sequence<1, 1, 0>>::type;

    if (verbose) {
      print_index_sequence(LayersUpToOrderThree{});
    }
    BOOST_CHECK(typeid(std::index_sequence<2, 2, 1>) ==
                typeid(LayersUpToOrderThree));
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
