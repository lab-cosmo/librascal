/**
 * @file test_sparse_kernels.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 Jan 2020
 *
 * @brief test the implementation of sparse similarity kernel classes
 *
 * @section LICENSE
 *
 * Copyright  2020 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_sparse_kernels.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {
  BOOST_AUTO_TEST_SUITE(sparse_kernels_test);

  using multiple_fixtures =
      boost::mpl::list<SparseKernelFixture<StrictNLSparseKernelFixture>>;

  /**
   * Tests if the compute functionality matches the size of atoms/structures
   * given as input.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_kernel_compute_test, Fix,
                                   multiple_fixtures, Fix) {
    using Calculator_t = typename Fix::Calculator_t;
    auto & kernels = Fix::kernels;
    auto & representations = Fix::representations;
    auto & collections = Fix::collections;

    // use all the features as sparse point
    std::vector<std::vector<std::vector<int>>> selected_ids;
    for (auto & collection : collections) {
      selected_ids.emplace_back();
      for (auto & manager : collection) {
        selected_ids.back().emplace_back();
        int ii{0};
        for (auto center : manager) {
          (void)center;
          selected_ids.back().back().push_back(ii);
          ++ii;
        }
      }
    }

    int i_collection{0};
    for (auto & collection : collections) {
      for (auto & representation : representations) {
        PseudoPointsBlockSparse<Calculator_t> sparse_points{};
        sparse_points.push_back(representation, collection,
                                selected_ids[i_collection]);
        for (auto & kernel : kernels) {
          auto mat = kernel.compute(representation, collection, sparse_points);

          if (Fix::verbose) {
            std::cout << "target_type=" << static_cast<int>(kernel.target_type)
                      << " mat.size=" << mat.size() << std::endl;
          }

          if (kernel.target_type == internal::TargetType::Structure) {
            BOOST_CHECK_EQUAL(mat.size(),
                              collection.size() * sparse_points.size());

          } else if (kernel.target_type == internal::TargetType::Atom) {
            int n_centers{0};
            for (auto & manager : collection) {
              n_centers += manager->size();
            }
            BOOST_CHECK_EQUAL(mat.size(), n_centers * sparse_points.size());
          }
        }
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
