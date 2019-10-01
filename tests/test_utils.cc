/**
 * file   test_utils.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2019
 *
 * @brief  test some utility functions
 *
 * Copyright  2019 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_utils.hh"

namespace rascal {
  BOOST_AUTO_TEST_SUITE(utils_test);

  /* ---------------------------------------------------------------------- */
  /**
   * Test the row norm sorting
   */
  BOOST_AUTO_TEST_CASE(cartesian_product_test) {
    bool verbose{false};
    std::vector<std::vector<int>> input{
      { { -1, 0 , 1 } , { 0 }, { -1 , 0 , 1 }}
    };
    std::vector<std::vector<int>> true_output{
      { { -1, 0 , -1 } , { -1, 0, 0 } , { -1, 0, 1  } , { 0, 0, -1 } ,
       { 0, 0, 0 } ,  { 0, 0, 1 } , { 1, 0, -1 } , { 1, 0, 0 } ,
       { 1, 0 ,1 }  }
    };

    auto output = internal::cartesian_product(input);
    BOOST_CHECK_EQUAL(true_output.size(), output.size());
    for (size_t iter{0}; iter < true_output.size(); iter++) {
      if (verbose) {
        std::cout << "true (";
        for (size_t ii{0}; ii < true_output[iter].size(); ii++) {
          std::cout << true_output[iter][ii] << ", ";
        }
        std::cout << ") "<< std::endl;

        std::cout << "test (";
        for (size_t ii{0}; ii < output[iter].size(); ii++) {
          std::cout << output[iter][ii] << ", ";
        }
        std::cout << ") "<< std::endl;

      }
      BOOST_CHECK_EQUAL_COLLECTIONS(
          true_output[iter].begin(),true_output[iter].end(),
          output[iter].begin(), output[iter].end());
    }

  }


  /* ---------------------------------------------------------------------- */
  /**
   * Test the row norm sorting
   */
  BOOST_AUTO_TEST_CASE(linear_index_test) {
    using internal::get_linear_index;
    using internal::get_multidim_index;
    bool verbose{false};
    std::array<int, 3> sizes{4,4,3};

    for (int ix{0}; ix < sizes[0]; ix++) {
      for (int iy{0}; iy < sizes[1]; iy++) {
        for (int iz{0}; iz < sizes[2]; iz++) {
          std::array<int, 3> input{ix,iy,iz};
          auto lin_idx = get_linear_index(sizes, input);
          auto output = get_multidim_index(sizes, lin_idx);
          if (verbose) {
            std::cout << "input: " << ix << ", "<<iy<< ", " << iz<< std::endl;
            std::cout << "output: " << output[0] << ", "
                      <<output[1] << ", " << output[2] << std::endl;
          }
          BOOST_CHECK_EQUAL_COLLECTIONS(
            input.begin(),input.end(),
            output.begin(), output.end());
        }
      }
    }
  }



  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
