/**
 * file   test_representation_manager.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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
  BOOST_AUTO_TEST_SUITE(representation_sorted_coulomb_test);

  /* ---------------------------------------------------------------------- */
  // test if it runs without seg fault and all
  using multiple_fixtures = boost::mpl::list<
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSortedCoulomb,
                          MultipleStructureSortedCoulomb,
                          Option::CMSortDistance>,
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSortedCoulomb,
                          MultipleStructureSortedCoulomb,
                          Option::CMSortRowNorm>>;
  // test if it reproduces the reference values
  using fixtures_ref_test = boost::mpl::list<
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSortedCoulomb,
                          SortedCoulombTestData,
                          Option::CMSortRowNorm>>;

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test,
                                   Fix, multiple_fixtures, Fix) {
    auto & managers = Fix::managers_strict;
    auto & representations = Fix::representations;
    auto & hypers = Fix::hypers;

    for (auto & manager : managers) {
      for (auto & hyper : hypers) {
        representations.emplace_back(manager, hyper);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_compute_test,
                                   Fix, multiple_fixtures, Fix) {
    auto & managers = Fix::managers_strict;
    auto & representations = Fix::representations;
    auto & hypers = Fix::hypers;
    for (auto& manager : managers) {
      for (auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
        representations.back().compute();
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_reference_test,
                                   Fix, fixtures_ref_test, Fix) {
    auto & managers = Fix::managers_strict;
    auto & representations = Fix::representations;
    auto & options = Fix::options;
    auto & hypers = Fix::hypers;
    auto & feature_matrices = Fix::feature_matrices;

    // Choose the data depending on the current options
    using Std2DArray_t = std::vector<std::vector<double>>;

    if (options[0] == Option::CMSortDistance) {
      hypers = Fix::data_sort_distance["hypers"];
      feature_matrices =
              Fix::data_sort_distance["feature_matrices"];
    } else if (options[0] == Option::CMSortRowNorm) {
      hypers = Fix::data_sort_rownorm["hypers"];
      feature_matrices =
              Fix::data_sort_rownorm["feature_matrices"];
    } else {
      hypers = Fix::data_sort_distance["hypers"];
      feature_matrices =
              Fix::data_sort_distance["feature_matrices"];
    }

    size_t manager_i{0};
    for (auto& manager : managers) {
      representations.emplace_back(manager,
                  hypers[manager_i].get<json>());
      representations.back().compute();
      const auto & ref_representation =
              feature_matrices[manager_i].get<Std2DArray_t>();

      const auto & test_representation =
                    representations.back().get_representation_full();

      if (manager_i == 2) {
      BOOST_CHECK_EQUAL(ref_representation.size(),
                        test_representation.rows());
      for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {
        BOOST_CHECK_EQUAL(ref_representation[row_i].size(),
                        test_representation.cols());
        for (size_t col_i{0}; col_i < ref_representation[row_i].size(); col_i++) {
          std::cout << std::abs(test_representation(row_i, col_i) - ref_representation[row_i][col_i]) << ", ";
          // BOOST_CHECK_EQUAL(ref_representation[row_i][col_i],
          //               test_representation(row_i, col_i));
        }
        std::cout << std::endl;
        if (row_i < 20) continue;
        if (row_i > 40) break;
      }
      std::cout << "#####################################################"<<std::endl;

      break;
      }

      manager_i += 1;
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
} // RASCAL
