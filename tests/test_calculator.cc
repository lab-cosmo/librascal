/**
 * file   test_calculator.cc
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
#include "test_calculator.hh"
#include "test_math.hh"  // for the gradient test

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

  using multiple_fixtures =
      boost::mpl::list<CalculatorFixture<MultipleStructureSortedCoulomb>,
                       CalculatorFixture<MultipleStructureSphericalExpansion>,
                       CalculatorFixture<MultipleStructureSphericalInvariants>,
                       CalculatorFixture<MultipleStructureSphericalCovariants>>;

  using fixtures_ref_test =
      boost::mpl::list<CalculatorFixture<SortedCoulombTestData>,
                       CalculatorFixture<SphericalExpansionTestData>,
                       CalculatorFixture<SphericalInvariantsTestData>,
                       CalculatorFixture<SphericalCovariantsTestData>>;

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the constructor runs and that the name is properly set
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test, Fix,
                                   multiple_fixtures, Fix) {
    auto & representations = Fix::representations;
    auto & representation_hypers = Fix::representation_hypers;
    bool verbose{false};

    for (auto & hyper : representation_hypers) {
      representations.emplace_back(hyper);
      auto & name{representations.back().get_name()};
      if (verbose) {
        std::cout << name << std::endl;
      }
    }
    // test the user defined name works
    for (auto & hyper : representation_hypers) {
      hyper["identifier"] = "my_representation";
      representations.emplace_back(hyper);
      auto & name{representations.back().get_name()};
      auto & prefix{representations.back().get_prefix()};
      BOOST_CHECK_EQUAL(prefix + std::string("my_representation"), name);
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test that the function get_dense_feature_matrix from managerCollection and
   * Property return the same dense matrix
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_dense_feature_comparison, Fix,
                                   multiple_fixtures, Fix) {
    using ManagerCollection_t =
        typename TypeHolderInjector<ManagerCollection,
                                    typename Fix::ManagerTypeList_t>::type;
    using Property_t = typename Fix::Property_t;

    bool verbose = false;

    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & representation_hypers = Fix::representation_hypers;
    int manager_i{0};
    for (auto & manager : managers) {
      for (auto & hyper : representation_hypers) {
        representations.emplace_back(hyper);
        representations.back().compute(manager);
        ManagerCollection_t collection{};
        auto & prop = manager->template get_validated_property_ref<Property_t>(
            representations.back().get_name());
        math::Matrix_t feat_prop = prop.get_dense_feature_matrix();
        collection.add_structure(manager);
        math::Matrix_t feat_col =
            collection.get_dense_feature_matrix(representations.back());

        BOOST_CHECK_EQUAL(feat_prop.rows(), feat_col.rows());
        BOOST_CHECK_EQUAL(feat_prop.cols(), feat_col.cols());
        for (int row_i{0}; row_i < feat_prop.rows(); row_i++) {
          for (int col_i{0}; col_i < feat_prop.cols(); ++col_i) {
            double diff =
                std::abs(feat_prop(row_i, col_i) - feat_col(row_i, col_i));

            BOOST_CHECK_LE(diff, 6e-12);
            if (verbose and diff > 6e-12) {
              std::cout << "manager_i=" << manager_i << " pos=" << row_i << ", "
                        << col_i << " \t " << feat_prop(row_i, col_i)
                        << "\t != " << feat_col(row_i, col_i) << std::endl;
            }
          }
        }
      }
      manager_i++;
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
    auto & representation_hypers = Fix::representation_hypers;
    for (auto & manager : managers) {
      for (auto & hyper : representation_hypers) {
        representations.emplace_back(hyper);
        representations.back().compute(manager);
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
    auto & verbose = Fix::verbose;
    using Property_t = typename Fix::Property_t;

    // Choose the data depending on the current options
    using Std2DArray_t = std::vector<std::vector<double>>;

    const auto & rep_infos{ref_data.at("rep_info").template get<json>()};
    // feature_matrices = data["feature_matrices"];

    size_t manager_i{0};
    for (auto & manager : managers) {
      for (const auto & rep_info : rep_infos.at(manager_i)) {
        const auto & representation_hypers =
            rep_info.at("hypers").template get<json>();
        const auto & ref_representation =
            rep_info.at("feature_matrix").template get<Std2DArray_t>();

        representations.emplace_back(representation_hypers);
        representations.back().compute(manager);
        auto property_name{representations.back().get_name()};
        auto && property{
            manager->template get_validated_property_ref<Property_t>(
                property_name)};
        auto test_representation = property.get_dense_feature_matrix();

        BOOST_CHECK_EQUAL(ref_representation.size(),
                          test_representation.rows());
        for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {
          BOOST_CHECK_EQUAL(ref_representation[row_i].size(),
                            test_representation.cols());
          for (size_t col_i{0}; col_i < ref_representation[row_i].size();
               ++col_i) {
            auto diff{std::abs(ref_representation[row_i][col_i] -
                               test_representation(row_i, col_i))};
            BOOST_CHECK_LE(diff, 6e-12);
            if (verbose and diff > 6e-12) {
              std::cout << "manager_i=" << manager_i << " pos=" << row_i << ", "
                        << col_i << " \t " << ref_representation[row_i][col_i]
                        << "\t != " << test_representation(row_i, col_i)
                        << std::endl;
            }
          }
        }
      }
      manager_i += 1;
    }
  }

  using fixtures_with_gradients =
      boost::mpl::list<CalculatorFixture<MultipleHypersSphericalExpansion>>;

  /**
   * Test the derivative of the GTO radial integral in the SphericalExpansion
   *
   * Doesn't depend much on the structuremanager or even the specific pair in
   * use; the nested loops below are just to pick out _a_ pair to supply as a
   * required argument to the radial integral functions.
   *
   * We do test a variety of hypers, though.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(spherical_expansion_radial_derivative, Fix,
                                   fixtures_with_gradients, Fix) {
    auto & managers = Fix::managers;
    auto & hypers = Fix::representation_hypers;
    // We need to explicitly specify a cluster ref type below - in this case,
    // it's for an atom pair (hence the 2)
    using ClusterRef_t = typename Fix::Manager_t::template ClusterRef<2>;
    using RadialIntegral_t =
        internal::RadialContribution<internal::RadialBasisType::GTO>;
    GradientTestFixture test_data{"reference_data/radial_derivative_test.json"};
    auto && it_manager{managers.front()->begin()};  // Need only one manager
    auto && atom{*it_manager};
    auto && it_atom{atom.begin()};
    auto && pair{*it_atom};  // Need only one (arbitrary) pair
    auto manager = managers.front();
    for (auto & hyper : hypers) {
      std::shared_ptr<RadialIntegral_t> radial_integral =
          std::make_shared<RadialIntegral_t>(hyper);
      // in C++17 the compiler would be able to deduce the template
      // arguments for itself >:/
      SphericalExpansionRadialDerivative<RadialIntegral_t, ClusterRef_t>
          calculator(radial_integral, pair);
      test_gradients(calculator, test_data);
    }
  }

  using simple_periodic_fixtures = boost::mpl::list<
      CalculatorFixture<SingleHypersSphericalRepresentation>,   // expension
      CalculatorFixture<SingleHypersSphericalRepresentation>>;  // invariants

  /**
   * Test the gradient of the SphericalExpansion representation on a few simple
   * crystal structures (single- and multi-species, primitive and supercells)
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(spherical_representation_gradients, Fix,
                                   simple_periodic_fixtures, Fix) {
    auto & managers = Fix::managers;
    auto & hypers = Fix::representation_hypers;
    auto & representations = Fix::representations;
    auto & structures = Fix::structures;
    auto filename_it = Fix::filenames.begin();

    for (auto manager : managers) {
      for (auto hyper : hypers) {
        hyper["compute_gradients"] = true;
        representations.emplace_back(hyper);
        structures.emplace_back();
        structures.back().set_structure(*filename_it);
        // The finite-difference tests don't work with periodic boundary
        // conditions -- moving one atom moves all its periodic images, too
        structures.back().pbc.setZero();
        RepresentationManagerGradientCalculator<typename Fix::Representation_t,
                                                typename Fix::Manager_t>
            calculator(representations.back(), manager, structures.back());
        RepresentationManagerGradientFixture<typename Fix::Representation_t,
                                             typename Fix::Manager_t>
            grad_fix("reference_data/spherical_expansion_gradient_test.json",
                     manager, calculator);
        if (grad_fix.verbosity >= GradientTestFixture::VerbosityValue::INFO) {
          std::cout << "Testing structure: " << *filename_it << std::endl;
        }
        do {
          test_gradients(grad_fix.get_calculator(), grad_fix);
          grad_fix.advance_center();
        } while (grad_fix.has_next());
      }
      ++filename_it;
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal