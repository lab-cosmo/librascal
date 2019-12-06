/**
 * @file   test_calculator.cc
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
#define BOOST_TEST_MODULE boost_test_sequence

#include "test_calculator.hh"

#include "test_math.hh"  // for the gradient test

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {
  /* ---------------------------------------------------------------------- */

  using multiple_fixtures =
      boost::mpl::list<CalculatorFixture<MultipleStructureSortedCoulomb<
                           MultipleStructureManagerNLStrictFixture>>,
                       CalculatorFixture<MultipleStructureSphericalExpansion<
                           MultipleStructureManagerNLCCStrictFixture>>,
                       CalculatorFixture<MultipleStructureSphericalInvariants<
                           MultipleStructureManagerNLCCStrictFixture>>,
                       CalculatorFixture<MultipleStructureSphericalCovariants<
                           MultipleStructureManagerNLCCStrictFixture>>>;

  using fixtures_ref_test =
      boost::mpl::list<CalculatorFixture<SortedCoulombTestData>,
                       CalculatorFixture<SphericalExpansionTestData>,
                       CalculatorFixture<SphericalInvariantsTestData>,
                       CalculatorFixture<SphericalCovariantsTestData>>;

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
   * Test that the function get_features from managerCollection and
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
        math::Matrix_t feat_prop = prop.get_features();
        collection.add_structure(manager);
        math::Matrix_t feat_col =
            collection.get_features(representations.back());

        BOOST_CHECK_EQUAL(feat_prop.rows(), feat_col.rows());
        BOOST_CHECK_EQUAL(feat_prop.cols(), feat_col.cols());
        double diff{0.};
        int size{0};
        for (int row_i{0}; row_i < feat_prop.rows(); row_i++) {
          for (int col_i{0}; col_i < feat_prop.cols(); ++col_i) {
            diff += std::abs(feat_prop(row_i, col_i) - feat_col(row_i, col_i));
            size += 1;

            if (verbose and diff / size > 6e-12) {
              std::cout << "manager_i=" << manager_i << " pos=" << row_i << ", "
                        << col_i << " \t " << feat_prop(row_i, col_i)
                        << "\t != " << feat_col(row_i, col_i) << std::endl;
            }
          }
        }
        diff /= size;
        BOOST_CHECK_LE(diff, 6e-12);
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
   * Test if the no center option takes out the centers
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_no_center_test, Fix,
                                   multiple_fixtures, Fix) {
    using ArrayB_t = typename AtomicStructure<3>::ArrayB_t;
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    using Property_t = typename Fix::Property_t;
    auto & hypers = Fix::representation_hypers;
    for (auto & manager : managers) {
      auto man = extract_underlying_manager<0>(manager);
      auto atomic_structure = man->get_atomic_structure();
      auto n_atoms = atomic_structure.get_number_of_atoms();
      atomic_structure.center_atoms_mask = ArrayB_t::Zero(n_atoms);
      auto i_atom1 = static_cast<int>(n_atoms / 2);
      atomic_structure.center_atoms_mask[i_atom1] = true;
      manager->update(atomic_structure);
      for (auto & hyper : hypers) {
        representations.emplace_back(hyper);
        std::string property_name{representations.back().get_name()};
        representations.back().compute(manager);
        auto prop{manager->template get_validated_property<Property_t>(
            property_name)};
        BOOST_CHECK_EQUAL(prop->get_nb_item(), 1);
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  using multiple_center_mask_fixtures = boost::mpl::list<
      // TODO(felix) For some reason the Sorted coulomb version is
      // susceptible to
      // differences in neighbour ordering so test will fail for a wrong
      // reason...
      // CalculatorFixture<MultipleStructureSortedCoulombCenterMask,
      //                       RepresentationManagerSortedCoulomb>,
      CalculatorFixture<MultipleStructureSphericalExpansion<
          MultipleStructureManagerNLCCStrictFixtureCenterMask>>,
      CalculatorFixture<MultipleStructureSphericalCovariants<
          MultipleStructureManagerNLCCStrictFixtureCenterMask>>,
      CalculatorFixture<MultipleStructureSphericalInvariants<
          MultipleStructureManagerNLCCStrictFixtureCenterMask>>>;

  /**
   * Test that selecting subsets of centers will give the same representation
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_center_mask_test, Fix,
                                   multiple_center_mask_fixtures, Fix) {
    bool verbose{false};
    auto & managers = Fix::managers;
    auto & hypers = Fix::representation_hypers;
    using Representation_t = typename Fix::Representation_t;
    using Property_t = typename Fix::Property_t;

    int n_manager{static_cast<int>(managers.size())};
    for (int i_manager{0}; i_manager < n_manager; i_manager += 2) {
      for (auto & hyper : hypers) {
        auto & manager = managers[i_manager];
        auto & manager_no_center = managers[i_manager + 1];
        auto center_atoms_mask =
            extract_underlying_manager<0>(manager_no_center)
                ->get_center_atoms_mask();

        if (verbose) {
          std::cout << "center_atoms_mask: " << center_atoms_mask.transpose()
                    << std::endl;
        }

        Representation_t representation{hyper};
        representation.compute(manager);
        representation.compute(manager_no_center);

        auto & prop = manager->template get_validated_property_ref<Property_t>(
            representation.get_name());
        math::Matrix_t rep_full = prop.get_features();

        auto & prop_no_center =
            manager_no_center->template get_validated_property_ref<Property_t>(
                representation.get_name());
        math::Matrix_t rep_no_center = prop_no_center.get_features();

        BOOST_CHECK_EQUAL(rep_full.cols(), rep_no_center.cols());
        BOOST_CHECK_EQUAL(center_atoms_mask.count(), rep_no_center.rows());

        if (verbose) {
          std::cout << "rep dim: " << rep_no_center.rows() << ", "
                    << rep_no_center.cols() << std::endl;
        }

        size_t i_no_center{0};
        for (size_t i_center{0}; i_center < manager_no_center->size();
             ++i_center) {
          if (center_atoms_mask(i_center)) {
            auto row_full = rep_full.row(i_center);
            auto row_no_center = rep_no_center.row(i_no_center);
            auto diff = (row_full - row_no_center).norm();
            BOOST_CHECK_LE(diff, math::DBL_FTOL);
            if (verbose) {
              std::cout << "Center idx: " << i_center << " Diff: " << diff
                        << std::endl;
            }
            i_no_center++;
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  using grad_sparse_fixtures =
      boost::mpl::list<CalculatorFixture<ComplexHypersSphericalInvariants>>;

  /**
   * Test if the keys for the PowerSpectrum's gradient are properly sparse on
   * the 1st center of CaCrP2O7_mvc-11955_symmetrized.json
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(sparse_grad_test, Fix, grad_sparse_fixtures,
                                   Fix) {
    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    using Representation_t = typename Fix::Representation_t;
    using PropGrad_t = typename Representation_t::template PropertyGradient_t<
        typename Fix::Manager_t>;
    auto & hypers = Fix::representation_hypers;
    const bool verbose{false};
    // assumes the structure is CaCrP2O7_mvc-11955_symmetrized.json
    // the keys of the center contribution (1st center)
    std::vector<std::vector<int>> all_keys{{8, 8},   {8, 15},  {8, 20},
                                           {15, 15}, {15, 20}, {20, 20}};
    // the list of keys for each neighbor (1st center)
    std::vector<std::vector<std::vector<int>>> neigh_keys{
        {{8, 15}, {15, 15}, {15, 20}}, {{8, 8}, {8, 15}, {8, 20}},
        {{8, 8}, {8, 15}, {8, 20}},    {{8, 8}, {8, 15}, {8, 20}},
        {{8, 8}, {8, 15}, {8, 20}},    {{8, 8}, {8, 15}, {8, 20}},
        {{8, 8}, {8, 15}, {8, 20}},    {{8, 8}, {8, 15}, {8, 20}},
        {{8, 8}, {8, 15}, {8, 20}},    {{8, 8}, {8, 15}, {8, 20}},
        {{8, 15}, {15, 15}, {15, 20}}};

    for (auto & manager : managers) {
      for (auto & hyper : hypers) {
        if (hyper["soap_type"] == "PowerSpectrum") {
          representations.emplace_back(hyper);
          std::string property_name{representations.back().get_gradient_name()};
          representations.back().compute(manager);
          auto && prop_grad{
              *manager->template get_validated_property<PropGrad_t>(
                  property_name)};
          for (auto center : manager) {
            auto ii_pair = center.get_atom_ii();
            auto keys_grad_center = prop_grad.get_keys(ii_pair);
            if (verbose) {
              std::cout << "Center gradient keys: ";
              for (auto key : keys_grad_center) {
                std::cout << "{";
                for (auto key_sp : key) {
                  std::cout << key_sp << ", ";
                }
                std::cout << "\b\b}, ";
              }
              std::cout << std::endl;
            }

            BOOST_TEST(keys_grad_center.size() == all_keys.size());
            for (size_t ii{0}; ii < keys_grad_center.size(); ii++) {
              BOOST_TEST(keys_grad_center[ii].size() == all_keys[ii].size());
              for (size_t jj{0}; jj < keys_grad_center[ii].size(); jj++) {
                BOOST_TEST(keys_grad_center[ii] == all_keys[ii],
                           boost::test_tools::per_element());
              }
            }
            int i_neigh{0};
            for (auto neigh : center.pairs()) {
              auto neigh_type = neigh.get_atom_type();
              auto keys_neigh = prop_grad[neigh].get_keys();
              if (verbose) {
                std::cout << "Neighbour " << neigh_type << " keys: ";
                for (auto key : keys_neigh) {
                  std::cout << "{";
                  for (auto key_sp : key) {
                    std::cout << key_sp << ", ";
                  }
                  std::cout << "\b\b}, ";
                }
                std::cout << std::endl;
              }

              BOOST_TEST(keys_neigh.size() == neigh_keys[i_neigh].size());
              for (size_t ii{0}; ii < keys_neigh.size(); ii++) {
                BOOST_TEST(keys_neigh[ii].size() == all_keys[ii].size());
                for (size_t jj{0}; jj < keys_neigh[ii].size(); jj++) {
                  BOOST_TEST(keys_neigh[ii] == neigh_keys[i_neigh][ii],
                             boost::test_tools::per_element());
                }
              }
              ++i_neigh;
            }
            break;  // only the first center
          }
        }
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
    auto verbose = true;
    using Property_t = typename Fix::Property_t;

    // Choose the data depending on the current options
    using Std2DArray_t = std::vector<std::vector<double>>;

    const auto & rep_infos{ref_data.at("rep_info").template get<json>()};

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
        auto test_representation = property.get_features();

        BOOST_CHECK_EQUAL(ref_representation.size(),
                          test_representation.rows());
        double avg_diff{0.};
        for (size_t row_i{0}; row_i < ref_representation.size(); row_i++) {
          BOOST_CHECK_EQUAL(ref_representation[row_i].size(),
                            test_representation.cols());
          for (size_t col_i{0}; col_i < ref_representation[row_i].size();
               ++col_i) {
            auto diff{std::abs(ref_representation[row_i][col_i] -
                               test_representation(row_i, col_i))};
            avg_diff += diff;
            if (verbose and diff > 6e-12) {
              std::cout << "manager_i=" << manager_i << " pos=" << row_i << ", "
                        << col_i << " \t " << ref_representation[row_i][col_i]
                        << "\t != " << test_representation(row_i, col_i)
                        << std::endl;
            }
          }
          BOOST_CHECK_LE(avg_diff / test_representation.size(), 6e-12);
        }
      }
      manager_i += 1;
    }
  }

  using fixtures_with_gradients = boost::mpl::list<
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::GTO,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::None>,
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::GTO,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::Interpolator>,
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::DVR,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::None>,
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::DVR,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::Interpolator>>;

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
    using RadialIntegral_t = typename Fix::RadialIntegral_t;

    GradientTestFixture test_data{"reference_data/tests_only/"
                                  "radial_derivative_test.json"};
    auto && it_manager{managers.front()->begin()};  // Need only one manager
    auto && atom{*it_manager};
    auto && it_atom{atom.pairs().begin()};
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

  using gradient_fixtures = boost::mpl::list<
      CalculatorFixture<
          SingleHypersSphericalExpansion<SimplePeriodicNLCCStrictFixture>>,
      CalculatorFixture<
          SingleHypersSphericalInvariants<SimplePeriodicNLCCStrictFixture>>>;

  /**
   * Test the gradient of the SphericalExpansion and SphericalInvariants
   * representation on a few simple crystal structures (single- and
   * multi-species, primitive and supercells)
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(spherical_representation_gradients, Fix,
                                   gradient_fixtures, Fix) {
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
        /* ---- grad-test-example-start1 ---- */
        RepresentationCalculatorGradientProvider<typename Fix::Representation_t,
                                                 typename Fix::Manager_t>
            provider(representations.back(), manager, structures.back());
        RepresentationCalculatorGradientFixture<typename Fix::Representation_t,
                                                typename Fix::Manager_t>
            grad_fix("reference_data/tests_only/"
                     "spherical_expansion_gradient_test.json",
                     manager, provider);
        /* ---- grad-test-example-end1 ---- */
        if (grad_fix.verbosity >= GradientTestFixture::VerbosityValue::INFO) {
          std::cout << "Testing structure: " << *filename_it << std::endl;
          std::cout << "With hypers: " << hyper << std::endl;
        }
        /* ---- grad-test-example-start2 ---- */
        do {
          test_gradients(grad_fix.get_provider(), grad_fix);
          grad_fix.advance_center();
        } while (grad_fix.has_next());
        /* ---- grad-test-example-end2 ---- */
      }
      ++filename_it;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Utility fixture used to compare representations computed with full and
   * half neighbor lists.
   */
  template <class CaculatorFixtureFull, class CaculatorFixtureHalf>
  struct MergeHalfAndFull : CaculatorFixtureFull, CaculatorFixtureHalf {
    using ParentFull = CaculatorFixtureFull;
    using ParentHalf = CaculatorFixtureHalf;
    using Representation_t = typename ParentFull::Representation_t;
    using Manager_t = typename ParentFull::Manager_t;
    using ManagerHalf_t = typename ParentHalf::Manager_t;
    using Prop_t = typename Representation_t::template Property_t<Manager_t>;
    using PropGrad_t =
        typename Representation_t::template PropertyGradient_t<Manager_t>;
    using PropHalf_t =
        typename Representation_t::template Property_t<ManagerHalf_t>;
    using PropGradHalf_t =
        typename Representation_t::template PropertyGradient_t<ManagerHalf_t>;

    MergeHalfAndFull() {
      for (auto hyper : ParentFull::representation_hypers) {
        hyper["compute_gradients"] = true;
        ParentFull::representations.emplace_back(hyper);
      }
    }

    ~MergeHalfAndFull() = default;
  };
  template <template <class> class RepresentationFixture>
  using RepFix_t =
      CalculatorFixture<RepresentationFixture<SimplePeriodicNLCCStrictFixture>>;
  template <template <class> class RepresentationFixture>
  using RepFixHalf_t = CalculatorFixture<
      RepresentationFixture<SimplePeriodicNLHalfCCStrictFixture>>;

  using gradient_half_fixtures = boost::mpl::list<
      MergeHalfAndFull<RepFix_t<SingleHypersSphericalExpansion>,
                       RepFixHalf_t<SingleHypersSphericalExpansion>>,
      MergeHalfAndFull<RepFix_t<SingleHypersSphericalInvariants>,
                       RepFixHalf_t<SingleHypersSphericalInvariants>>>;

  /**
   * Test the representation gradients computed with a half neighbor list
   * against the full neighbor list implementation.
   * There are some discrepancies between the descriptors, their gradients and
   * their half NL conterparts arising from numerical errors most likely.
   *
   * Extensive tests on a subset of the structures showed that large relative
   * errors can arise when both the reference and the test values are small
   * (~1e-18) so we use a threshold (epsilon ~ 1e-15) so that if the absolute
   * value of reference and test values are within epsilon the relative error
   * is effectively 0. This is justified because these coefficient will be
   * either taken as is or multiplied together (-> power 2, 3, 4...) which
   * will squash these values even more.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(half_full_representation_test, Fix,
                                   gradient_half_fixtures, Fix) {
    using Prop_t = typename Fix::Prop_t;
    using PropHalf_t = typename Fix::PropHalf_t;
    using PropGrad_t = typename Fix::PropGrad_t;
    using PropGradHalf_t = typename Fix::PropGradHalf_t;
    auto & managers = Fix::ParentFull::managers;
    auto & managers_half = Fix::ParentHalf::managers;
    auto & representations = Fix::ParentFull::representations;
    const bool verbose{false};
    // relative error threshold
    const double delta{4e-7};
    // range of zero
    const double epsilon{1e-15};

    // both manager should refer to the same structures
    BOOST_TEST(managers.size() == managers_half.size());
    if (managers.size() != managers_half.size()) {
      throw std::runtime_error("managers.size() != managers_half.size()");
    }
    for (size_t i_manager{0}; i_manager < managers.size(); ++i_manager) {
      for (auto & representation : representations) {
        auto & manager = managers[i_manager];
        auto & manager_half = managers_half[i_manager];
        representation.compute(manager);
        representation.compute(manager_half);

        auto && rep_vectors{*manager->template get_property_ptr<Prop_t>(
            representation.get_name())};
        auto && rep_vectors_half{
            *manager_half->template get_property_ptr<PropHalf_t>(
                representation.get_name())};

        auto && rep_vector_gradients{
            *manager->template get_property_ptr<PropGrad_t>(
                representation.get_gradient_name())};
        auto && rep_vector_gradients_half{
            *manager_half->template get_property_ptr<PropGradHalf_t>(
                representation.get_gradient_name())};

        size_t center_count{0};
        for (auto center : manager) {
          // compare the representation coefficients
          auto diff_rep_m{math::relative_error(
              rep_vectors.get_dense_row(center),
              rep_vectors_half.get_dense_row(center), delta, epsilon)};
          double diff_rep = diff_rep_m.maxCoeff();
          BOOST_TEST(diff_rep < delta);
          if (verbose and diff_rep > delta) {
            std::cout << "========================= rep" << std::endl;
            std::cout << "Center " << center.get_index();
            std::cout << " of type " << center.get_atom_type()
                      << " max rel diff: " << diff_rep << std::endl;
            std::cout << "Full: " << std::endl
                      << rep_vectors.get_dense_row(center).transpose();
            std::cout << std::endl;
            std::cout << "Half: " << std::endl
                      << rep_vectors_half.get_dense_row(center).transpose();
            std::cout << std::endl;
          }

          auto half_it = manager_half->get_iterator_at(center_count, 0);
          auto half_center = *(half_it);
          auto ii_pair = center.get_atom_ii();
          auto ii_half_pair = half_center.get_atom_ii();

          // compare the representation gradient coefficients at the ii pair
          auto diff_rep_grad_center_m{math::relative_error(
              rep_vector_gradients.get_dense_row(ii_pair),
              rep_vector_gradients_half.get_dense_row(ii_half_pair), delta,
              epsilon)};
          double diff_rep_grad_center = diff_rep_grad_center_m.maxCoeff();
          BOOST_TEST(diff_rep_grad_center < delta);
          if (verbose and diff_rep_grad_center > delta) {
            std::cout << "================ rep_grad_center" << std::endl;
            std::cout << "Center " << center.get_index();
            std::cout << " of type " << center.get_atom_type()
                      << " max rel diff: " << diff_rep_grad_center << std::endl;
            std::cout
                << "Full: " << std::endl
                << rep_vector_gradients.get_dense_row(ii_pair).transpose();
            std::cout << std::endl;
            std::cout << "Half: " << std::endl
                      << rep_vector_gradients_half.get_dense_row(ii_half_pair)
                             .transpose();
            std::cout << std::endl;
          }
          size_t neigh_count{0};
          for (auto neigh : center.pairs()) {
            auto neigh_type = neigh.get_atom_type();
            auto tags = neigh.get_atom_tag_list();
            if (tags[1] <= tags[0]) {
              continue;
            }

            auto half_neigh_it = half_center.pairs().begin();
            for (size_t ii{0}; ii < neigh_count; ii++) {
              ++half_neigh_it;
            }
            auto half_neigh = *(half_neigh_it);
            // compare the representation gradient coefficients at ij pair
            auto diff_rep_grad_neigh_m{math::relative_error(
                rep_vector_gradients.get_dense_row(neigh),
                rep_vector_gradients_half.get_dense_row(half_neigh), delta,
                epsilon)};
            double diff_rep_grad_neigh = diff_rep_grad_neigh_m.maxCoeff();
            BOOST_TEST(diff_rep_grad_neigh < delta);
            if (verbose and diff_rep_grad_neigh > delta) {
              std::cout << "================== rep_grad_neigh" << std::endl;
              std::cout << "Center " << center.get_index();
              std::cout << " of type " << center.get_atom_type() << std::endl;
              std::cout << "Neighbour " << neigh_type << " tags: ";
              std::cout << "(";
              for (auto tag : tags) {
                std::cout << tag << ", ";
              }
              std::cout << "\b\b) "
                        << "max rel diff: " << diff_rep_grad_neigh << std::endl;
              std::cout
                  << "Full: " << std::endl
                  << rep_vector_gradients.get_dense_row(neigh).transpose();
              std::cout << std::endl;
              std::cout << "Half: " << std::endl
                        << rep_vector_gradients_half.get_dense_row(half_neigh)
                               .transpose();
              std::cout << std::endl;
            }
            neigh_count++;
          }  // neigh
          center_count++;
        }  // center
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
