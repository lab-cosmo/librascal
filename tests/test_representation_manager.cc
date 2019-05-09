/**
 * file   test_representation_manager.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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
#include "test_representation_manager.hh"

namespace rascal {
  BOOST_AUTO_TEST_SUITE(representation_test);

  /* ---------------------------------------------------------------------- */
  /**
   * Test the row norm sorting
   */
  BOOST_AUTO_TEST_CASE(rownorm_sort_test) {
    Eigen::MatrixXd test_matrix(4, 5);
    // clang-format off
    test_matrix << 0, 6, 1, 4, 3,
                   0, 7, 2, 5, 4,
                   1, 8, 3, 6, 2,
                   2, 9, 4, 7, 1;
    // clang-format on
    Eigen::MatrixXd true_order(5, 1);
    // use of stable sort so 2 goes before 4
    true_order << 0, 1, 3, 2, 4;

    auto test_order =
        internal::SortCoulomMatrix<internal::CMSortAlgorithm::RowNorm>::
            get_coulomb_matrix_sorting_order(test_matrix, test_matrix);

    for (auto idx_i{0}; idx_i < true_order.size(); ++idx_i) {
      BOOST_CHECK_EQUAL(true_order(idx_i), test_order[idx_i].first);
    }
  }

  /**
   * Test the distance from the central atom sorting.
   * assumes the center is on row 0.
   */
  BOOST_AUTO_TEST_CASE(distance_sort_test) {
    Eigen::MatrixXd test_matrix(4, 4);
    // clang-format off
    test_matrix << 0.        , 1.68624958, 1.43774399, 1.12522187,
                   1.68624958,         0.,  1.6850887, 1.15322292,
                   1.43774399,  1.6850887,         0., 0.98009938,
                   1.12522187, 1.15322292, 0.98009938,         0.;
    // clang-format on
    Eigen::MatrixXd true_order(4, 1);
    // use of stable sort so 2 goes before 4
    true_order << 0, 3, 2, 1;

    auto test_order =
        internal::SortCoulomMatrix<internal::CMSortAlgorithm::Distance>::
            get_coulomb_matrix_sorting_order(test_matrix, test_matrix);

    for (auto idx_i{0}; idx_i < true_order.size(); ++idx_i) {
      BOOST_CHECK_EQUAL(true_order(idx_i), test_order[idx_i].first);
    }
  }

  /* ---------------------------------------------------------------------- */

  using multiple_fixtures = boost::mpl::list<
      RepresentationFixture<MultipleStructureSortedCoulomb,
                            RepresentationManagerSortedCoulomb>,
      RepresentationFixture<MultipleStructureSphericalExpansion,
                            RepresentationManagerSphericalExpansion>,
      RepresentationFixture<MultipleStructureSOAP, RepresentationManagerSOAP>>;

  using fixtures_ref_test = boost::mpl::list<RepresentationFixture<
      SortedCoulombTestData, RepresentationManagerSortedCoulomb>>;

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the constructor runs
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test, Fix,
                                   multiple_fixtures, Fix) {
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & hypers = Fix::hypers;

    for (auto & manager : managers) {
      for (auto & hyper : hypers) {
        representations.emplace_back(manager, hyper);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the compute function runs
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_compute_test, Fix,
                                   multiple_fixtures, Fix) {
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & hypers = Fix::hypers;
    for (auto & manager : managers) {
      for (auto & hyper : hypers) {
        representations.emplace_back(manager, hyper);
        representations.back().compute();
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the representation computed is equal to a reference from a file
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_reference_test, Fix,
                                   fixtures_ref_test, Fix) {
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & ref_data = Fix::ref_data;

    // Choose the data depending on the current options
    using Std2DArray_t = std::vector<std::vector<double>>;

    const auto & rep_infos{ref_data.at("rep_info").template get<json>()};
    // feature_matrices = data["feature_matrices"];

    size_t manager_i{0};
    for (auto & manager : managers) {
      for (const auto & rep_info : rep_infos.at(manager_i)) {
        const auto & hypers = rep_info.at("hypers").template get<json>();
        const auto & ref_representation =
            rep_info.at("feature_matrix").template get<Std2DArray_t>();

        representations.emplace_back(manager, hypers);
        representations.back().compute();
        auto aa{hypers.dump(2)};
        const auto & test_representation =
            representations.back().get_representation_full();

        BOOST_CHECK_EQUAL(ref_representation.size(),
                          test_representation.rows());
        for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {
          BOOST_CHECK_EQUAL(ref_representation[row_i].size(),
                            test_representation.cols());

          for (size_t col_i{0}; col_i < ref_representation[row_i].size();
               ++col_i) {
            auto diff{std::abs(ref_representation[row_i][col_i] -
                               test_representation(row_i, col_i))};
            BOOST_CHECK_LE(diff, 1e-12);
          }
        }
      }
      manager_i += 1;
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

  /* ---------------------------------------------------------------------- */

  /* Tests specific to the spherical expansion representation
   * Test the SOAP representation
   */
  BOOST_AUTO_TEST_SUITE(representation_blocksparse_specific_tests);

  using fixtures_ref_test = boost::mpl::list<
      RepresentationFixture<SphericalExpansionTestData,
                            RepresentationManagerSphericalExpansion>,
      RepresentationFixture<SOAPTestData, RepresentationManagerSOAP>>;

  /*
   * Test if the representation computed is equal to a reference from a file
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_reference_test, Fix,
                                   fixtures_ref_test, Fix) {
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & ref_data = Fix::ref_data;

    // Choose the data depending on the current options
    // using Std2DArray_t = std::vector<std::vector<double>>;
    using Std2DArray_t = std::vector<std::vector<double>>;

    const auto & data{ref_data.at("rep_info").template get<json>()};
    // feature_matrices = data["feature_matrices"];

    size_t manager_i{0};
    for (auto & manager : managers) {
      for (const auto & config : data.at(manager_i)) {
        const auto & hypers = config.at("hypers").template get<json>();
        const auto & ref_representation =
            config.at("feature_matrix").template get<Std2DArray_t>();

        representations.emplace_back(manager, hypers);
        representations.back().compute();

        // TODO(felix) quick fix of something that will disappear soon
        FeatureManagerBlockSparse<double> features{
            representations.back().get_feature_size(), hypers};
        features.push_back(representations.back());
        auto test_representation = features.get_feature_matrix_dense();

        auto n_feature{test_representation.rows()};
        auto n_center{test_representation.cols()};
        BOOST_CHECK_EQUAL(ref_representation.size(), n_feature);
        for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {
          BOOST_CHECK_EQUAL(ref_representation[row_i].size(), n_center);
          for (size_t col_i{0}; col_i < ref_representation[row_i].size();
               ++col_i) {
            auto diff{std::abs(ref_representation[row_i][col_i] -
                               test_representation(row_i, col_i))};
            BOOST_CHECK_LE(diff, 1e-12);
          }
        }
      }
      manager_i += 1;
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
