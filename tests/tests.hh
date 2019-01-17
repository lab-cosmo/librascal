#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#ifndef TESTS_TESTS_HH_
#define TESTS_TESTS_HH_

#include <iostream>
#include <iomanip>
#include <boost/test/unit_test.hpp>
//#define BOOST_TEST_DYN_LINK

namespace rascal {

  const double tol = 1e-14 * 100;  // it's in percent

}  // namespace rascal

#endif  // TESTS_TESTS_HH_
