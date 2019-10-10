#include "test_cluster_ref_key.hh"


namespace rascal {
  BOOST_AUTO_TEST_SUITE(tests_cluster_ref_key);

  /**
   * tests the get_min_layer function
   */
  BOOST_AUTO_TEST_CASE( test_get_min_layer ) {
    // LayersByOrder using only one root manager
    const size_t active_layer_root_manager{get_min_layer<1>(typename std::index_sequence<0>{})};    
    BOOST_CHECK_EQUAL(0, active_layer_root_manager);
    // LayersByOrder for `AdaptorFullNeighbourList`
    const size_t active_layer_order_one_adaptor{get_min_layer<1>(typename std::index_sequence<1,0>{})};
    BOOST_CHECK_EQUAL(1, active_layer_order_one_adaptor);
    // LayersByOrder for stacking an `AdaptorMaxOrder` on an `AdaptorStrict`
    const size_t active_layer_increase_max_order{get_min_layer<3>(typename std::index_sequence<1,1,0>{})};
    BOOST_CHECK_EQUAL(0, active_layer_increase_max_order);
  }

  /**
   * tests the get_layer function
   */
  BOOST_FIXTURE_TEST_CASE( test_get_layer, LayerFixture) {
    BOOST_CHECK(get_layer<1>(order_three_seq{}) == 1);
    BOOST_CHECK(get_layer<2>(order_three_seq{}) == 1);
    BOOST_CHECK(get_layer<3>(order_three_seq{}) == 0);
  }

  /**
   * tests the get_layers function
   */
  BOOST_FIXTURE_TEST_CASE( test_get_layers, LayerFixture) {
    std::array<size_t, 1> arr_order_one = get_layers(order_one_seq{});
    std::array<size_t, 1> arr_order_one_ref = {1};
    BOOST_CHECK(arr_order_one == arr_order_one_ref);

    std::array<size_t, 2> arr_order_two = get_layers(order_two_seq{});
    std::array<size_t, 2> arr_order_two_ref = {1,1};
    BOOST_CHECK(arr_order_two == arr_order_two_ref);

    std::array<size_t, 3> arr_order_three = get_layers(order_three_seq{});
    std::array<size_t, 3> arr_order_three_ref = {1,1,0};
    BOOST_CHECK(arr_order_three == arr_order_three_ref);
    print_index_sequence(std::make_index_sequence<5>{});
  }

  /**
   * tests the LayerExtender struct 
   */
  BOOST_AUTO_TEST_CASE( test_layer_extender ) {
    bool verbose{true};
    
    using LayersUpToOrderTwo = typename LayerExtender<std::index_sequence<1>>::type;
    using LayersUpToOrderThree = typename LayerExtender<LayersUpToOrderTwo>::type;

    if (verbose) {
      print_index_sequence(LayersUpToOrderTwo{});
    }
    BOOST_CHECK( typeid(std::index_sequence<1, 0>) == typeid(LayersUpToOrderTwo) );

    if (verbose) {
      print_index_sequence(LayersUpToOrderThree{});
    }
    BOOST_CHECK( typeid(std::index_sequence<1, 0, 0>) == typeid(LayersUpToOrderThree) );
  }

  /**
   * tests the LayerIncreaser struct 
   */
  BOOST_AUTO_TEST_CASE( test_layer_increaser ) {
    bool verbose{true};
    using LayersUpToOrderThree = typename LayerIncreaser<std::index_sequence<1,1,0>>::type;

    if (verbose) {
      print_index_sequence(LayersUpToOrderThree{});
    }
    BOOST_CHECK( typeid(std::index_sequence<2, 2, 1>) == typeid(LayersUpToOrderThree) );
  }
  
  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
