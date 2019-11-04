/**
 * @file   test_utils.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2019
 *
 * @brief Test some utility functions
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

// Detects which compiler is used
#if defined(__clang__)
#define CLANG_COMPILER
#elif ((defined(__GNUC__) || defined(__GNUG__)) && (__GNUC__ > 6))  // NOLINT
#define GCC_COMPILER_7_AND_UPPER
#elif ((defined(__GNUC__) || defined(__GNUG__)) && (__GNUC__ <= 6))  // NOLINT
#define GCC_COMPILER_6_AND_LOWER
#endif

#include "test_utils.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(test_utils);
  /* ---------------------------------------------------------------------- */
  /**
   * Test that the internal::GetTypeNameHelper function returns what we expect
   */
  BOOST_AUTO_TEST_CASE(type_demangling_test) {
    bool verbose{false};

    // test template with numbers
    std::string ref_clang{"ClusterRefKey_1_2"};
    std::string ref_gcc_6_and_lower{"ClusterRefKey_1ul_2ul"};
    std::string ref_gcc_7_and_upper{"ClusterRefKey_1_2"};
    std::string test{internal::GetTypeName<ClusterRefKey<1, 2>>()};
#if defined(GCC_COMPILER_7_AND_UPPER)
    BOOST_CHECK_EQUAL(ref_gcc_7_and_upper, test);
#elif defined(GCC_COMPILER_6_AND_LOWER)
    BOOST_CHECK_EQUAL(ref_gcc_6_and_lower, test);
#elif defined(CLANG_COMPILER)
    BOOST_CHECK_EQUAL(ref_clang, test);
#endif
    if (verbose) {
      std::cout << "ref_clang:  " << ref_clang << std::endl;
      std::cout << "ref_gcc_6_and_lower:  " << ref_gcc_6_and_lower << std::endl;
      std::cout << "ref_gcc_7_and_upper:  " << ref_gcc_7_and_upper << std::endl;
      std::cout << "test: " << test << std::endl;
      std::cout << "pretty_function: "
                << internal::PrettyFunction<ClusterRefKey<1, 2>>::get()
                << std::endl;
    }

    // test template with class
    std::string ref_simple{"StructureManager_StructureManagerCenters"};
    std::string test_simple{
        internal::GetTypeName<StructureManager<StructureManagerCenters>>()};
    BOOST_CHECK_EQUAL(ref_simple, test_simple);
    if (verbose) {
      std::cout << "ref_simple:  " << ref_simple << std::endl;
      std::cout << "test_simple:  " << test_simple << std::endl;
      std::cout << "pretty_function: "
                << internal::PrettyFunction<
                       StructureManager<StructureManagerCenters>>::get()
                << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
