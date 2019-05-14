/**
 * file    test_feature_manager.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief  test  managers
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

#include "tests.hh"
#include "test_feature_manager.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(feature_dense_test);
  /* ---------------------------------------------------------------------- */
  // TODO(felix) define more test that could be streamlined
  // gets a list of fixtures for all the different possible structure managers
  using multiple_fixtures =
      boost::mpl::list<FeatureFixture<double, FeatureManagerDense,
                                      RepresentationManagerSortedCoulomb,
                                      MultipleStructureSortedCoulomb>,
                       FeatureFixture<float, FeatureManagerDense,
                                      RepresentationManagerSortedCoulomb,
                                      MultipleStructureSortedCoulomb>>;

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the Fixture with multiple structures builds
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_setup_test, Fix, multiple_fixtures,
                                   Fix) {}

  /* ---------------------------------------------------------------------- */
  /**
   * Test the construction of the feature manager
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test, Fix,
                                   multiple_fixtures, Fix) {
    auto & features = Fix::features;
    auto & hypers = Fix::hypers;
    auto & n_feature = Fix::n_feature;
    for (auto & hyper : hypers) {
      features.emplace_back(n_feature, hyper);
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test pushing back the feature matrices of several representations in the
   * feature manager and check if the data is properly set.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(representation_aggregate_test, Fix,
                                   multiple_fixtures, Fix) {
    using Precision_t = typename Fix::Precision_t;
    auto & features = Fix::features;
    auto & hypers = Fix::hypers;
    auto & n_feature = Fix::n_feature;
    auto & n_center = Fix::n_center;

    // build the feature managers. only the 1st
    // one will be used
    for (auto & hyper : hypers) {
      features.emplace_back(n_feature, hyper);
      features.back().reserve(n_center);
    }

    auto & representations = Fix::representations;

    std::vector<std::vector<Precision_t>> original_data{};
    std::vector<int> n_centers{};
    std::vector<int> n_features{};
    // extract the feature matrices in a ref vector
    for (auto & representation : representations) {
      auto && raw_data{representation.get_representation_raw_data()};
      original_data.emplace_back(raw_data.begin(), raw_data.end());
      n_centers.push_back(representation.get_center_size());
      n_features.push_back(representation.get_feature_size());
    }

    // move the features into the feature manager
    for (auto & representation : representations) {
      features.front().push_back(representation);
    }

    // check if the features have been moved properly
    auto feature_matrix = features.front().get_feature_matrix();
    int stride{0};
    for (size_t it{0}; it < original_data.size(); ++it) {
      for (int icenter{0}; icenter < n_centers[it]; ++icenter) {
        for (int ifeature{0}; ifeature < n_features[it]; ++ifeature) {
          int lin_idx{icenter * n_features[it] + ifeature};
          double diff{original_data[it][lin_idx] -
                      feature_matrix(ifeature, stride + icenter)};
          BOOST_CHECK_LE(diff, 1e-14);
        }
      }
      stride += n_centers[it];
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

  BOOST_AUTO_TEST_SUITE(feature_block_sparse_test);

  using multiple_fixtures_sparse = boost::mpl::list<
      SparseFeatureFixture<double, FeatureManagerBlockSparse,
                           RepresentationManagerSphericalExpansion,
                           MultipleStructureSphericalExpansion>,
      SparseFeatureFixture<float, FeatureManagerBlockSparse,
                           RepresentationManagerSphericalExpansion,
                           MultipleStructureSphericalExpansion>>;

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the Fixture with multiple structures builds
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_setup_test, Fix,
                                   multiple_fixtures_sparse, Fix) {}
  /* ---------------------------------------------------------------------- */
  /**
   * Test the construction of the feature manager
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test, Fix,
                                   multiple_fixtures_sparse, Fix) {
    auto & features = Fix::features;
    auto & hypers = Fix::hypers;
    auto & inner_sizes = Fix::inner_sizes;
    int ii{0};
    for (auto & hyper : hypers) {
      features.emplace_back(inner_sizes[ii], hyper);
      ii++;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test pushing back the feature matrices of several representations in the
   * feature manager and check if the data is properly set.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(representation_aggregate_test, Fix,
                                   multiple_fixtures_sparse, Fix) {
    bool verbose{false};
    auto & features = Fix::features;
    auto & hypers = Fix::hypers;
    auto & inner_sizes = Fix::inner_sizes;
    auto & representations = Fix::representations;

    // build the feature managers. only the 1st
    // one will be used
    for (size_t i_hyper{0}; i_hyper < hypers.size(); i_hyper++) {
      features.emplace_back(inner_sizes[i_hyper], hypers[i_hyper]);
    }

    using Data_t = typename Fix::Representation_t::SparseProperty_t::Data_t;
    std::vector<Data_t> original_data{};
    // extract the feature matrices in a ref vector
    for (size_t i_hyper{0}; i_hyper < hypers.size(); i_hyper++) {
      original_data.emplace_back();
      for (auto & representation : representations[i_hyper]) {
        auto & raw_data{representation.get_representation_sparse_raw_data()};
        original_data.back().insert(original_data.back().end(),
                                    raw_data.begin(), raw_data.end());
      }
    }

    // move the features into the feature manager
    for (size_t i_hyper{0}; i_hyper < hypers.size(); i_hyper++) {
      for (size_t i_rep{0}; i_rep < representations[i_hyper].size(); i_rep++) {
        features[i_hyper].push_back(representations[i_hyper][i_rep]);
      }
    }

    // check if the features have been moved properly
    for (size_t i_hyper{0}; i_hyper < hypers.size(); i_hyper++) {
      auto feature_matrix = features[i_hyper].get_feature_matrix_dense();
      auto & inner_size{inner_sizes[i_hyper]};
      std::set<std::vector<int>> unique_keys{};
      for (size_t i_center{0}; i_center < original_data[i_hyper].size();
           ++i_center) {
        for (const auto & element : original_data[i_hyper][i_center]) {
          unique_keys.emplace(element.first);
        }
      }
      for (size_t i_center{0}; i_center < original_data[i_hyper].size();
           ++i_center) {
        auto && datas{original_data[i_hyper][i_center]};
        int i_count{0};
        double diff{0};
        int N{0};
        if (verbose)
          std::cout << "Center: " << i_center << std::endl;
        for (auto & key : unique_keys) {
          // if the key exist
          if (datas.count(key) == 1) {
            if (verbose)
              std::cout << "Key: " << key[0] << std::endl;
            auto && data = datas[key];
            for (int i_col{0}; i_col < data.cols(); i_col++) {
              for (int i_row{0}; i_row < data.rows(); i_row++) {
                diff += data(i_row, i_col) - feature_matrix(i_count, i_center);
                N += 1;
                if (verbose) {
                  if (diff > 1e-12) {
                    std::cout << diff << ", ";
                    std::cout << data(i_row, i_col) << ", ";
                    std::cout << feature_matrix(i_count, i_center) << "/ ";
                  }
                }

                i_count++;
              }
              if (verbose)
                std::cout << std::endl;
            }
          } else {
            i_count += inner_size;
          }
        }
        BOOST_CHECK_LE(diff / N, 1e-10);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
