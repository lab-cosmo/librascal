/**
 * file    test_feature_manager.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief  test  managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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
  using multiple_fixtures = boost::mpl::list<
    FeatureFixture<double, FeatureManagerDense,
                   StructureManagerCenters,
                   RepresentationManagerSortedCoulomb,
                   TestFeatureData, Option::CMSortDistance>,
    FeatureFixture<double, FeatureManagerDense,
                   StructureManagerCenters,
                   RepresentationManagerSortedCoulomb,
                   TestFeatureData, Option::CMSortRowNorm>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_setup_test,
                                   Fix, multiple_fixtures, Fix) {
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test,
                                   Fix, multiple_fixtures, Fix) {
    auto& features = Fix::features;
    auto& hypers = Fix::hypers;
    auto& Nfeature = Fix::Nfeature;
    for (auto& hyper : hypers) {
      features.emplace_back(Nfeature, hyper);
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(representation_aggregate_test,
                                   Fix, multiple_fixtures, Fix) {
    auto& features = Fix::features;
    auto& hypers = Fix::hypers;
    auto& Nfeature = Fix::Nfeature;
    auto& Ncenter = Fix::Ncenter;
    for (auto& hyper : hypers) {
      features.emplace_back(Nfeature, hyper);
      features.back().reserve(Ncenter);
    }

    auto& representations = Fix::representations;
    for (auto& representation : representations) {
      features.front().push_back(representation);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
} // namespace rascal
