

#  if defined(__INTEL_COMPILER)
//#    pragma warning ( disable : 383 )
#  elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Weffc++"
#  elif (defined(__GNUC__) || defined(__GNUG__))
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Weffc++"
#  endif
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#  if defined(__INTEL_COMPILER)
//#    pragma warning ( disable : 383 )
#  elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#    pragma clang diagnostic pop
#    pragma clang diagnostic ignored "-Weffc++"
#  elif (defined(__GNUC__) || defined(__GNUG__))
#    pragma GCC diagnostic pop
#    pragma GCC diagnostic ignored "-Weffc++"
#  endif

#ifndef TESTS_H
#define TESTS_H

#include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

namespace proteus {

  const double tol = 1e-14*100; //it's in percent

}  // proteus


#endif /* TESTS_H */
