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


#include "test_utils.hh"

namespace rascal {


  BOOST_AUTO_TEST_SUITE(test_utils);
  /* ---------------------------------------------------------------------- */
  /**
   * Test that the internal::GetTypeNameHelper function
   */
  BOOST_AUTO_TEST_CASE(type_demangling_test) {
    bool verbose{false};
    std::string ref{"rascal::ClusterRefKey<1, 2>"};
    std::string test{internal::GetTypeNameHelper<ClusterRefKey<1,2>>::GetTypeName()};
    BOOST_CHECK_EQUAL(ref, test);
    if (verbose) {
      std::cout << "ref:  "<< ref << std::endl;
      std::cout << "test: "<<test << std::endl;
      std::cout << "pretty_function: " << internal::PrettyFunction<ClusterRefKey<1,2>>::get() << std::endl;
    }
  }


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
