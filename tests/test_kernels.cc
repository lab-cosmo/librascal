/**
 * @file test_kernels.cc
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
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "test_kernels.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {
  BOOST_AUTO_TEST_SUITE(kernels_test);

  using multiple_ref_fixtures =
      boost::mpl::list<KernelFixture<DataSphericalInvariantsKernelFixture>>;

  /**
   * Tests if the wrong target_type is catched correctly.
   */
  BOOST_FIXTURE_TEST_CASE(kernel_target_type_test,
                          KernelFixture<DataSphericalInvariantsKernelFixture>) {
    auto kernel_hyper = this->ParentA::kernel_hypers.at(0);
    kernel_hyper["target_type"] = "this_target_type_does_not_exist";
    BOOST_CHECK_THROW(auto kernel_wrong_target_type = Kernel(kernel_hyper),
                      std::runtime_error);
    kernel_hyper["target_type"] = "structure";
    BOOST_CHECK_THROW(auto kernel_wrong_target_type = Kernel(kernel_hyper),
                      std::runtime_error);
    kernel_hyper.erase("target_type");
    BOOST_CHECK_THROW(auto kernel_no_target_type = Kernel(kernel_hyper),
                      std::runtime_error);
  }

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
    const double delta{1e-10};
    BOOST_CHECK_EQUAL(collections.size(), ref_data.size());
    for (size_t i_collection{0}; i_collection < ref_data.size();
         ++i_collection) {
      auto & collection = collections[i_collection];
      BOOST_CHECK_EQUAL(kernels.size(), ref_data[i_collection].size());
      for (size_t i_rep{0}; i_rep < ref_data[i_collection].size(); ++i_rep) {
        auto ref_mat = ref_data[i_collection][i_rep]["kernel_matrix"]
                           .template get<math::Matrix_t>();
        auto & rep = representations[i_rep];
        auto & kernel = kernels[i_rep];
        rep.compute(collection);
        auto mat = kernel.compute(rep, collection, collection);
        auto diff_m{math::relative_error(ref_mat, mat, delta)};
        double diff = diff_m.maxCoeff();
        BOOST_TEST(diff < delta);
        if (diff > delta) {
          std::cout << ref_data[i_collection][i_rep]["hypers_kernel"].dump()
                    << std::endl;
          std::cout << ref_data[i_collection][i_rep]["hypers_rep"].dump()
                    << std::endl;
          std::cout << mat.row(0) << std::endl;
          std::cout << ref_mat.row(0) << std::endl;
        }

        auto mat_sym = kernel.compute(rep, collection);
        diff_m = math::relative_error(ref_mat, mat_sym, delta);
        diff = diff_m.maxCoeff();
        BOOST_TEST(diff < delta);
      }
    }
  }

  /**
   * Tests that the Spherical invariant give the same kernel with different
   * expansion_by_species_method
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(species_method_test, Fix,
                                   multiple_ref_fixtures, Fix) {
    using Calculator_t = typename Fix::Calculator_t;
    auto & kernels = Fix::kernels;
    auto & representation_hypers = Fix::ParentB::representation_hypers;
    auto & collections = Fix::collections;
    auto & ref_data = Fix::ParentA::ref_data;
    const double delta{1e-10};
    BOOST_CHECK_EQUAL(collections.size(), ref_data.size());
    for (size_t i_collection{0}; i_collection < ref_data.size();
         ++i_collection) {
      auto & collection = collections[i_collection];
      BOOST_CHECK_EQUAL(kernels.size(), ref_data[i_collection].size());
      for (size_t i_rep{0}; i_rep < ref_data[i_collection].size(); ++i_rep) {
        auto ref_mat = ref_data[i_collection][i_rep]["kernel_matrix"]
                           .template get<math::Matrix_t>();

        auto & rep_hypers = representation_hypers[i_rep];
        auto & kernel = kernels[i_rep];
        rep_hypers["expansion_by_species_method"] = "environment wise";
        Calculator_t rep_env_wise{rep_hypers};
        rep_env_wise.compute(collection);
        auto mat = kernel.compute(rep_env_wise, collection);
        auto diff_m{math::relative_error(ref_mat, mat, delta)};
        double diff = diff_m.maxCoeff();
        BOOST_TEST(diff < delta);

        rep_hypers["expansion_by_species_method"] = "structure wise";
        Calculator_t rep_structure_wise{rep_hypers};
        rep_structure_wise.compute(collection);
        mat = kernel.compute(rep_structure_wise, collection);
        diff_m = math::relative_error(ref_mat, mat, delta);
        diff = diff_m.maxCoeff();
        BOOST_TEST(diff < delta);

        rep_hypers["expansion_by_species_method"] = "user defined";
        std::vector<int> global_species{1, 6, 7, 8};
        rep_hypers["global_species"] = global_species;
        Calculator_t rep_user_defined{rep_hypers};
        rep_user_defined.compute(collection);
        mat = kernel.compute(rep_user_defined, collection);
        diff_m = math::relative_error(ref_mat, mat, delta);
        diff = diff_m.maxCoeff();
        BOOST_TEST(diff < delta);
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
