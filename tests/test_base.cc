
#include "tests.hh"
#include "module.hh"

#include <boost/mpl/list.hpp>


namespace proteus {

  BOOST_AUTO_TEST_SUITE(base_tests);

  BOOST_AUTO_TEST_CASE(base_test) {
    BOOST_CHECK_EQUAL(1, 2-1);
  }

  BOOST_AUTO_TEST_CASE(f_test) {
    BOOST_TEST(f(12) == 2, "Should have been 2");
  }


  template <int Dim>
  struct DemoTestFixture {

    static constexpr int dim(){return Dim;}

    DemoTestFixture()
      :val{Dim}
    {}


    int val;
  };

  using fixtures = boost::mpl::list<DemoTestFixture<2>, DemoTestFixture<3>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(templated_basic_fixture_test, Fix, fixtures, Fix) {
    BOOST_TEST(Fix::val == Fix::dim());
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // proteus
