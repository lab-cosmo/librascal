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
        SparsePointsBlockSparse<Calculator_t> sparse_points{};
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

  /* ---------------------------------------------------------------------- */
  /**
   * Utility fixture used to compare representations with sparsification
   */

  struct SparseKernelGradFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    using Manager_t = typename ManagerTypeHolder_t::type;
    using ManagerCollection_t =
        typename TypeHolderInjector<ManagerCollection,
                                    ManagerTypeHolder_t::type_list>::type;
    using Structure_t = AtomicStructure<3>;
    using Representation_t = CalculatorSphericalInvariants;
    using SparsePoints_t = SparsePointsBlockSparse<Representation_t>;
    using Kernel_t = SparseKernel;
    using Prop_t = typename Representation_t::template Property_t<Manager_t>;
    using PropGrad_t =
        typename Representation_t::template PropertyGradient_t<Manager_t>;
    json inputs{};

    SparseKernelGradFixture() {
      this->inputs =
          json_io::load("reference_data/tests_only/sparse_kernel_inputs.json");
    }

    ~SparseKernelGradFixture() = default;
  };

  using sparse_grad_fixtures = boost::mpl::list<SparseKernelGradFixture>;

  /**
   * Test the analytical kernel gradients against numerical kernel gradients.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(grad_test, Fix, sparse_grad_fixtures, Fix) {
    // using Manager_t = typename Fix::Manager_t;
    using ManagerCollection_t = typename Fix::ManagerCollection_t;
    using Representation_t = typename Fix::Representation_t;
    using Kernel_t = typename Fix::Kernel_t;
    using SparsePoints_t = typename Fix::SparsePoints_t;

    const auto & inputs = Fix::inputs;

    const bool verbose{true};
    // relative error threshold
    const double delta{5e-6};
    // range of zero
    const double epsilon{1e-15};

    for (const auto & input : inputs) {
      // extract inputs
      std::string filename{input.at("filename").template get<std::string>()};
      json adaptors_input = input.at("adaptors").template get<json>();
      json calculator_input = input.at("calculator").template get<json>();
      json kernel_input = input.at("kernel").template get<json>();
      auto selected_ids = input.at("selected_ids")
                              .template get<std::vector<std::vector<int>>>();
      // initialize classes
      Kernel_t kernel{kernel_input};
      kernel_input.at("target_type") = "Atom";
      Kernel_t kernel_num{kernel_input};
      ManagerCollection_t managers{adaptors_input};
      SparsePoints_t sparse_points{};
      Representation_t representation{calculator_input};
      // load structures, compute representation and fill sparse points
      managers.add_structures(filename, 0,
                              input.at("n_structures").template get<int>());
      representation.compute(managers);
      sparse_points.push_back(representation, managers, selected_ids);
      // compute kernel gradients
      auto KNM_der{
          kernel.compute_derivative(representation, managers, sparse_points)};
      auto KNM_num_der{compute_numerical_kernel_gradients(
          kernel_num, representation, managers, sparse_points,
          input.at("h").template get<double>())};
      auto diff = math::relative_error(KNM_der, KNM_num_der, delta, epsilon);
      int col_max{0}, row_max{0};
      double max_rel_diff{diff.maxCoeff(&row_max, &col_max)};
      BOOST_TEST(max_rel_diff < delta);
      if (verbose and max_rel_diff > delta) {
        std::cout << filename << std::endl;
        std::cout << adaptors_input.dump() << std::endl;
        std::cout << calculator_input.dump() << std::endl;
        std::cout << kernel_input.dump() << std::endl;
        std::cout << diff.row(row_max) << std::endl;
        std::cout << "============================" << std::endl;
        std::cout << KNM_der.row(row_max) << std::endl;
        std::cout << "============================" << std::endl;
        std::cout << KNM_num_der.row(row_max) << std::endl;
        std::cout << "============================" << std::endl;
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
