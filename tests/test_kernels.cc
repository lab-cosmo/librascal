/**
 * file test_kernels.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 June 2019
 *
 * @brief test the implementation of similarity kernel classes
 *
 * @section LICENSE
 *
 * Copyright  2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_kernels.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {
  BOOST_AUTO_TEST_SUITE(kernels_test);

  using multiple_ref_fixtures =
      boost::mpl::list<KernelFixture<DataSphericalInvariantsKernelFixture>>;

  /**
   * Tests if the compute functionality agrees with the results of the reference
   * data.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(kernel_ref_data_test, Fix,
                                   multiple_ref_fixtures, Fix) {
    auto & kernels = Fix::kernels;
    auto & representations = Fix::representations;
    auto & collections = Fix::collections;
    auto & ref_data = Fix::ParentA::ref_data;

    using Std2DArray_t = std::vector<std::vector<double>>;

    BOOST_CHECK_EQUAL(collections.size(), ref_data.size());
    for (size_t i_collection{0}; i_collection < ref_data.size();
         ++i_collection) {
      auto & collection = collections[i_collection];
      BOOST_CHECK_EQUAL(kernels.size(), ref_data[i_collection].size());
      for (size_t i_rep{0}; i_rep < ref_data[i_collection].size(); ++i_rep) {
        auto & rep = representations[i_rep];
        auto & kernel = kernels[i_rep];
        rep.compute(collection);
        auto mat = kernel.compute(rep, collection, collection);
        auto ref_mat = ref_data[i_collection][i_rep]["kernel_matrix"]
                           .template get<Std2DArray_t>();
        double error{0.};
        for (int i_row{0}; i_row < mat.rows(); ++i_row) {
          for (int i_col{0}; i_col < mat.cols(); ++i_col) {
            auto error_ = std::abs(mat(i_row, i_col) - ref_mat[i_row][i_col]);
            if (error < error_) {
              error = error_;
            }
          }
        }
        BOOST_CHECK_LE(error, 6e-13);
      }
    }
  }

  using multiple_fixtures =
      boost::mpl::list<KernelFixture<StrictNLKernelFixture>>;

  /**
   * Tests if the compute functionality matches the size of atoms/structures
   * given as input.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_kernel_compute_test, Fix,
                                   multiple_fixtures, Fix) {
    auto & kernels = Fix::kernels;
    auto & representations = Fix::representations;
    auto & collections = Fix::collections;

    for (auto & collection : collections) {
      for (auto & representation : representations) {
        for (auto & kernel : kernels) {
          auto mat = kernel.compute(representation, collection, collection);

          if (Fix::verbose) {
            std::cout << "target_type=" << static_cast<int>(kernel.target_type)
                      << " mat.size=" << mat.size() << std::endl;
          }

          if (kernel.target_type == internal::TargetType::Structure) {
            BOOST_CHECK_EQUAL(mat.size(),
                              collection.size() * collection.size());

          } else if (kernel.target_type == internal::TargetType::Atom) {
            int n_centers{0};
            for (auto & manager : collection) {
              n_centers += manager->size();
            }
            BOOST_CHECK_EQUAL(mat.size(), n_centers * n_centers);
          }
        }
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
