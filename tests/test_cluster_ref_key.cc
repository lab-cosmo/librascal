#include "tests.hh"
#include "structure_managers/cluster_ref_key.hh"

#include <tuple>

namespace rascal {
  BOOST_AUTO_TEST_SUITE(cluster_ref_key);

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_CASE( test_case1 )
  {
    const size_t cluster_layer{compute_cluster_layer<1>(typename std::index_sequence<0>{})};
    BOOST_CHECK_EQUAL(0, cluster_layer);
  }
  BOOST_AUTO_TEST_CASE( test_case2 )
  {
    const size_t cluster_layer{compute_cluster_layer<2>(typename std::index_sequence<0,1>{})};
    BOOST_CHECK_EQUAL(1, cluster_layer);
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
