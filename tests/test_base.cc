#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>


namespace proteus {

  BOOST_AUTO_TEST_CASE(base_test) {
    BOOST_CHECK_EQUAL(1, 2-1);
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
    BOOST_CHECK_EQUAL(Fix::val, Fix::dim());
  }

}  // proteus
