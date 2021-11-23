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

  using spherical_fixtures =
      boost::mpl::list<CalculatorFixture<MultipleStructureSphericalExpansion<
                           MultipleStructureManagerNLCCStrictFixture>>,
                       CalculatorFixture<MultipleStructureSphericalInvariants<
                           MultipleStructureManagerNLCCStrictFixture>>>;

  using fixtures_ref_test =
      boost::mpl::list<CalculatorFixture<SortedCoulombTestData>,
                       CalculatorFixture<SphericalExpansionTestData>,
                       CalculatorFixture<SphericalInvariantsTestData>,
                       CalculatorFixture<SphericalCovariantsTestData>>;

  BOOST_AUTO_TEST_SUITE(representation_test);

  // the interaction cutoff is stored under different names for different
  // representations
  double extract_interaction_cutoff_from_representation_hyper(json hyper) {
    double representation_cutoff{0};
    if (hyper.find("cutoff_function") != hyper.end()) {
      // spherical representations
      representation_cutoff =
          hyper.at("cutoff_function").at("cutoff").at("value");
    } else if (hyper.find("interaction_cutoff") != hyper.end()) {
      // coulomb representation
      representation_cutoff = hyper.at("interaction_cutoff");
    } else {
      std::stringstream err_str{};
      err_str << "interaction cutoff of not found in representation "
              << "hyperparameters. "
              << "Please check representation hyperparameters: " << std::endl
              << hyper << std::endl;
      throw std::runtime_error(err_str.str());
    }
    return representation_cutoff;
  }

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

    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & representation_hypers = Fix::representation_hypers;
    int manager_i{0};
    std::set<int> species_list;
    for (auto & manager : managers) {
      for (auto center : manager) {
        species_list.insert(center.get_atom_type());
      }

      for (auto & hyper : representation_hypers) {
        double representation_cutoff{
            extract_interaction_cutoff_from_representation_hyper(hyper)};
        if (manager->get_cutoff() == representation_cutoff) {
          representations.emplace_back(hyper);
          representations.back().compute(manager);
          ManagerCollection_t collection{};
          auto & prop = *manager->template get_property<Property_t>(
              representations.back().get_name(), true);
          math::Matrix_t feat_prop = prop.get_features();
          collection.add_structure(manager);
          math::Matrix_t feat_col =
              collection.get_features(representations.back());

          BOOST_CHECK_EQUAL(feat_prop.rows(), feat_col.rows());
          BOOST_CHECK_EQUAL(feat_prop.cols(), feat_col.cols());
          auto diff_rep{math::relative_error(feat_prop, feat_col)};
          BOOST_CHECK_LE(diff_rep.maxCoeff(), 6e-12);
        }
      }
      manager_i++;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(spherical_number_feature_comparison, Fix,
                                   spherical_fixtures, Fix) {
    using Property_t = typename Fix::Property_t;

    auto & managers = Fix::managers;
    auto & representations = Fix::representations;
    auto & representation_hypers = Fix::representation_hypers;
    int manager_i{0};
    std::set<int> species_list;
    for (auto & manager : managers) {
      species_list.clear();
      for (auto center : manager) {
        species_list.insert(center.get_atom_type());
      }
      for (auto & hyper : representation_hypers) {
        double representation_cutoff{
            extract_interaction_cutoff_from_representation_hyper(hyper)};
        if (manager->get_cutoff() == representation_cutoff) {
          representations.emplace_back(hyper);
          representations.back().compute(manager);
          auto & prop = *manager->template get_property<Property_t>(
              representations.back().get_name(), true);
          math::Matrix_t feat_prop = prop.get_features();
          BOOST_CHECK_EQUAL(
              feat_prop.cols(),
              representations.back().get_num_coefficients(species_list.size()));
        }
      }
      manager_i++;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the compute function runs and
   * Test (de)serialization using the nlohmann::json format
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(serialization_test, Fix, multiple_fixtures,
                                   Fix) {
    using Representation_t = typename Fix::Representation_t;
    using Property_t = typename Fix::Property_t;
    auto & managers = Fix::managers;
    auto & representation_hypers = Fix::representation_hypers;

    for (auto & hyper : representation_hypers) {
      Representation_t rep{hyper};
      // serialize the calculator
      json j = rep;
      // deserialize the calculator
      Representation_t rep_test{j.template get<Representation_t>()};
      bool testing{rep == rep_test};
      // test serialization
      BOOST_TEST(testing == true);

      for (auto & manager : managers) {
        // does the compute function runs without errors ?
        rep.compute(manager);
        auto prop_ref =
            manager->template get_property<Property_t>(rep.get_name(), true);
        auto feat_ref = prop_ref->get_features();

        rep_test.compute(manager);
        auto prop_test = manager->template get_property<Property_t>(
            rep_test.get_name(), true);
        auto feat_test = prop_test->get_features();

        auto diff = math::relative_error(feat_ref, feat_test);
        // are the 2 identical calculators giving the same features ?
        BOOST_TEST(diff.maxCoeff() < 1e-12);
        break;  // a bit too long otherwise
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
        double representation_cutoff{
            extract_interaction_cutoff_from_representation_hyper(hyper)};
        if (manager->get_cutoff() == representation_cutoff) {
          representations.emplace_back(hyper);
          std::string property_name{representations.back().get_name()};
          representations.back().compute(manager);
          auto prop{manager->template get_property<Property_t>(
              representations.back().get_name(), true)};
          BOOST_CHECK_EQUAL(prop->get_nb_item(), 1);
        }
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
        double representation_cutoff{
            extract_interaction_cutoff_from_representation_hyper(hyper)};
        if (manager->get_cutoff() == representation_cutoff) {
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

          auto & prop = *manager->template get_property<Property_t>(
              representation.get_name(), true);
          math::Matrix_t rep_full = prop.get_features();

          auto & prop_no_center =
              *manager_no_center->template get_property<Property_t>(
                  representation.get_name(), true);
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
              *manager->template get_property<PropGrad_t>(property_name, true)};
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
              BOOST_TEST(keys_grad_center[ii] == all_keys[ii],
                         boost::test_tools::per_element());
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
              if (not hyper["normalize"]) {
                BOOST_TEST(keys_neigh.size() == neigh_keys[i_neigh].size());
                for (size_t ii{0}; ii < keys_neigh.size(); ii++) {
                  BOOST_TEST(keys_neigh[ii] == neigh_keys[i_neigh][ii],
                             boost::test_tools::per_element());
                }
              } else {
                BOOST_TEST(keys_neigh.size() == all_keys.size());
                for (size_t ii{0}; ii < keys_neigh.size(); ii++) {
                  BOOST_TEST(keys_neigh[ii] == all_keys[ii],
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
    const double delta{3e-6};
    const double epsilon{1e-14};
    using Property_t = typename Fix::Property_t;

    const auto & rep_infos{ref_data.at("rep_info").template get<json>()};

    size_t manager_i{0};
    for (auto & manager : managers) {
      for (const auto & rep_info : rep_infos.at(manager_i)) {
        const auto & representation_hypers =
            rep_info.at("hypers").template get<json>();
        const auto & ref_representation =
            rep_info.at("feature_matrix").template get<math::Matrix_t>();

        representations.emplace_back(representation_hypers);
        representations.back().compute(manager);
        auto property_name{representations.back().get_name()};
        auto && property{
            *manager->template get_property<Property_t>(property_name, true)};
        auto test_representation = property.get_features();

        auto diff_rep{math::relative_error(
            ref_representation, test_representation, delta, epsilon)};
        if (diff_rep.maxCoeff() > delta) {
          std::cout << "representation_hypers " << representation_hypers
                    << std::endl;
          BOOST_CHECK_LE(diff_rep.maxCoeff(), delta);
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
                                   internal::OptimizationType::Spline>,
      RadialIntegralHandlerFixture<
          MultipleHypersSphericalExpansion, internal::RadialBasisType::GTO,
          internal::AtomicSmearingType::Constant,
          internal::OptimizationType::RadialDimReductionSpline>,
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::DVR,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::None>,
      RadialIntegralHandlerFixture<MultipleHypersSphericalExpansion,
                                   internal::RadialBasisType::DVR,
                                   internal::AtomicSmearingType::Constant,
                                   internal::OptimizationType::Spline>,
      RadialIntegralHandlerFixture<
          MultipleHypersSphericalExpansion, internal::RadialBasisType::DVR,
          internal::AtomicSmearingType::Constant,
          internal::OptimizationType::RadialDimReductionSpline>>;

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
        // for peformance reasons we do the SphericalInvariants test only for
        // the first hyper, since it builds on top of the SphericalExpansion and
        // the SphericalExpansion is tested for all hypers, we do not lose
        // coverage
        if (hyper.count("soap_type")) {
          break;
        }
      }
      ++filename_it;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Utility fixture used to compare representations with sparsification
   */

  struct SparsificationSphericalInvariantsFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    using Structure_t = AtomicStructure<3>;
    using Representation_t = CalculatorSphericalInvariants;

    json factory_args{};
    std::vector<Structure_t> structures{};
    std::vector<json> representation_hypers{};
    std::vector<json> representation_sparse_hypers{};
    std::vector<std::vector<int>> selected_features_global_ids{};

    SparsificationSphericalInvariantsFixture() {
      json datas =
          json_io::load("reference_data/tests_only/sparsification_inputs.json");
      for (const auto & data : datas) {
        json structure{{"filename", data.at("filename").get<std::string>()}};
        auto adaptors = data.at("hypers").at("adaptors").get<json>();
        auto rep_hyp = data.at("hypers").at("rep").get<json>();
        auto rep_sparse_hyp = data.at("hypers").at("rep_sparse").get<json>();
        auto selected_features_global_ids =
            rep_sparse_hyp.at("coefficient_subselection")
                .at("selected_features_global_ids")
                .get<std::vector<int>>();
        json parameters;
        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;
        this->factory_args.emplace_back(parameters);
        this->representation_hypers.emplace_back(rep_hyp);
        this->representation_sparse_hypers.emplace_back(rep_sparse_hyp);
        this->selected_features_global_ids.emplace_back(
            selected_features_global_ids);
      }
    }

    ~SparsificationSphericalInvariantsFixture() = default;
  };

  using sparsification_fixtures = boost::mpl::list<
      CalculatorFixture<SparsificationSphericalInvariantsFixture>>;

  /**
   * Test the sparsification of the representation
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(sparsification_test, Fix,
                                   sparsification_fixtures, Fix) {
    using Manager_t = typename Fix::Manager_t;
    using Representation_t = typename Fix::Representation_t;
    using Prop_t = typename Representation_t::template Property_t<Manager_t>;
    using PropGrad_t =
        typename Representation_t::template PropertyGradient_t<Manager_t>;

    auto & managers = Fix::managers;
    auto & representation_hypers = Fix::representation_hypers;
    auto & representation_sparse_hypers = Fix::representation_sparse_hypers;
    auto & selected_features_global_ids = Fix::selected_features_global_ids;

    const bool verbose{true};
    // relative error threshold
    const double delta{4e-7};
    // range of zero
    const double epsilon{1e-15};
    // to iterate through gradients info
    int grandients_info_i_row{0};

    for (size_t i_manager{0}; i_manager < managers.size(); ++i_manager) {
      auto manager = managers[i_manager];
      Representation_t representation{representation_hypers[i_manager]};
      representation.compute(manager);

      Representation_t representation_sparse{
          representation_sparse_hypers[i_manager]};
      representation_sparse.compute(manager);

      auto & selected_ids = selected_features_global_ids[i_manager];

      auto & rep{
          *manager->template get_property<Prop_t>(representation.get_name())};
      auto & rep_sparse{*manager->template get_property<Prop_t>(
          representation_sparse.get_name())};

      // test that the underlying data of sparse version is a subset of the
      // original coefficient
      for (auto center : manager) {
        auto rep_r = rep[center].get_full_vector();
        auto rep_sparse_r = rep_sparse[center].get_full_vector();
        for (size_t i_feat{0}; i_feat < selected_ids.size(); ++i_feat) {
          double rel_err{math::relative_error(rep_r[selected_ids[i_feat]],
                                              rep_sparse_r[i_feat], delta,
                                              epsilon)};
          BOOST_TEST(rel_err < delta);
          if (verbose and rel_err > delta) {
            std::cout << "############## Man: " << i_manager
                      << " Center tag: " << center.get_atom_tag() << std::endl;
            std::cout << "REF:  " << rep_r.transpose() << std::endl;
            std::cout << "TEST: " << rep_sparse_r.transpose() << std::endl;
          }
        }
      }

      // test that get_features properly output the underlying data
      auto X = rep.get_features();
      auto X_sparse = rep_sparse.get_features();
      for (size_t i_feat{0}; i_feat < selected_ids.size(); ++i_feat) {
        auto rel_err_m{math::relative_error(
            X.col(selected_ids[i_feat]), X_sparse.col(i_feat), delta, epsilon)};
        double rel_err = rel_err_m.maxCoeff();
        BOOST_TEST(rel_err < delta);
      }

      // test that the sparsified gradients are a subset of the normal gradients
      auto & rep_grad{*manager->template get_property<PropGrad_t>(
          representation.get_gradient_name())};
      auto & rep_grad_sparse{*manager->template get_property<PropGrad_t>(
          representation_sparse.get_gradient_name())};
      auto X_grad = rep_grad.get_features_gradient();
      auto X_grad_sparse = rep_grad_sparse.get_features_gradient();
      for (size_t i_feat{0}; i_feat < selected_ids.size(); ++i_feat) {
        auto rel_err_m{math::relative_error(X_grad.col(selected_ids[i_feat]),
                                            X_grad_sparse.col(i_feat), delta,
                                            epsilon)};
        double rel_err = rel_err_m.maxCoeff();
        BOOST_TEST(rel_err < delta);
        if (verbose and rel_err > delta) {
          std::cout << "############## Man: " << i_manager << " col: " << i_feat
                    << std::endl;
          std::cout << "REF:  " << X_grad.col(selected_ids[i_feat]).transpose()
                    << std::endl;
          std::cout << "TEST: " << X_grad_sparse.col(i_feat).transpose()
                    << std::endl;
        }

        Eigen::Matrix<int, Eigen::Dynamic, 5> gradients_info{
            manager->get_gradients_info()};
        BOOST_TEST(gradients_info.rows() * 3 == X_grad.rows(),
                   "REF X_grad number of rows:  "
                       << X_grad.rows()
                       << " not equal to TEST gradients_info number of rows*3: "
                       << gradients_info.rows() * 3);
        grandients_info_i_row = 0;
        // equal check for the individual manager
        for (auto center : manager) {
          for (auto pair : center.pairs_with_self_pair()) {
            // because we do not create gradients info over manager collection
            // i_frame is 0
            BOOST_TEST(gradients_info(grandients_info_i_row, 0) == 0);
            BOOST_TEST(gradients_info(grandients_info_i_row, 1) ==
                       center.get_atom_tag());
            BOOST_TEST(gradients_info(grandients_info_i_row, 2) ==
                       pair.get_atom_j().get_atom_tag());
            BOOST_TEST(gradients_info(grandients_info_i_row, 3) ==
                       center.get_atom_type());
            BOOST_TEST(gradients_info(grandients_info_i_row, 4) ==
                       pair.get_atom_type());
            grandients_info_i_row++;
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Utility fixture used to compare representations computed with full and
   * half neighbor lists.
   */

  /** Contains some simple periodic structures for testing complicated things
   *  like gradients
   */
  struct SimpleFullFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    using Structure_t = AtomicStructure<3>;

    SimpleFullFixture() {
      for (size_t i_filename{0}; i_filename < filenames.size(); ++i_filename) {
        json parameters;
        json structure{{"filename", filenames[i_filename]}};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                 {"initialization_arguments",
                  {{"cutoff", cutoffs[i_filename]}, {"skin", cutoff_skin}}}};
        json ad1b{{"name", "AdaptorCenterContribution"},
                  {"initialization_arguments", {}}};
        json ad2{
            {"name", "AdaptorStrict"},
            {"initialization_arguments", {{"cutoff", cutoffs[i_filename]}}}};
        adaptors.emplace_back(ad1);
        adaptors.emplace_back(ad1b);
        adaptors.emplace_back(ad2);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);

        this->representation_hypers.emplace_back();
        auto fc_hyper = json::object(
            {{"type", "ShiftedCosine"},
             {"cutoff", {{"value", cutoffs[i_filename]}, {"unit", "AA"}}},
             {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}});

        for (auto & ri_hyp : this->radial_contribution_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyper;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.back().push_back(rep_hyp);
            }
          }
        }
      }
    }

    ~SimpleFullFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/inputs/diamond_2atom.json",
        "reference_data/inputs/diamond_2atom_distorted.json",
        "reference_data/inputs/SiCGe_wurtzite_like.json",
        "reference_data/inputs/methane.json"};
    const std::vector<double> cutoffs{{1.2, 1.2, 1.5, 4.}};
    const double cutoff_skin{0.};

    json factory_args{};
    std::vector<Structure_t> structures{};

    std::vector<std::vector<json>> representation_hypers{};
    std::vector<json> fc_hypers{};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{{{"max_radial", 2},
                                  {"max_angular", 2},
                                  {"normalize", true},
                                  {"soap_type", "PowerSpectrum"},
                                  {"compute_gradients", true}},
                                 {{"max_radial", 3},
                                  {"max_angular", 0},
                                  {"normalize", true},
                                  {"soap_type", "RadialSpectrum"},
                                  {"compute_gradients", true}}};
  };

  /** Contains some simple periodic structures for testing complicated things
   *  like gradients
   *  Should match the
   */
  struct SimpleHalfFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorHalfList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    using Structure_t = AtomicStructure<3>;

    SimpleHalfFixture() {
      for (size_t i_filename{0}; i_filename < filenames.size(); ++i_filename) {
        json parameters;
        json structure{{"filename", filenames[i_filename]}};
        json adaptors;
        json ad1a{
            {"name", "AdaptorNeighbourList"},
            {"initialization_arguments", {{"cutoff", cutoffs[i_filename]}}}};
        json ad1b{{"name", "AdaptorHalfList"},
                  {"initialization_arguments", {}}};
        json ad1c{{"name", "AdaptorCenterContribution"},
                  {"initialization_arguments", {}}};
        json ad2{
            {"name", "AdaptorStrict"},
            {"initialization_arguments", {{"cutoff", cutoffs[i_filename]}}}};

        adaptors.emplace_back(ad1a);
        adaptors.emplace_back(ad1b);
        adaptors.emplace_back(ad1c);
        adaptors.emplace_back(ad2);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);
      }
    }

    ~SimpleHalfFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/inputs/diamond_2atom.json",
        "reference_data/inputs/diamond_2atom_distorted.json",
        "reference_data/inputs/SiCGe_wurtzite_like.json",
        "reference_data/inputs/methane.json"};
    const std::vector<double> cutoffs{{1.2, 1.2, 1.5, 4.}};
    const double cutoff_skin{0.};

    json factory_args{};
    std::vector<Structure_t> structures{};
  };

  template <class FixtureFull, class FixtureHalf, class Calculator>
  struct MergeHalfAndFull : MultipleStructureFixture<FixtureFull>,
                            MultipleStructureFixture<FixtureHalf> {
    using ParentFull = MultipleStructureFixture<FixtureFull>;
    using ParentHalf = MultipleStructureFixture<FixtureHalf>;
    using Representation_t = Calculator;
    using Manager_t = typename ParentFull::Manager_t;
    using ManagerHalf_t = typename ParentHalf::Manager_t;
    using Prop_t = typename Representation_t::template Property_t<Manager_t>;
    using PropGrad_t =
        typename Representation_t::template PropertyGradient_t<Manager_t>;
    using PropHalf_t =
        typename Representation_t::template Property_t<ManagerHalf_t>;
    using PropGradHalf_t =
        typename Representation_t::template PropertyGradient_t<ManagerHalf_t>;

    MergeHalfAndFull() {}

    ~MergeHalfAndFull() = default;
  };

  using gradient_half_fixtures =
      boost::mpl::list<MergeHalfAndFull<SimpleFullFixture, SimpleHalfFixture,
                                        CalculatorSphericalExpansion>,
                       MergeHalfAndFull<SimpleFullFixture, SimpleHalfFixture,
                                        CalculatorSphericalInvariants>>;
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
    using Representation_t = typename Fix::Representation_t;
    auto & managers = Fix::ParentFull::managers;
    auto & managers_half = Fix::ParentHalf::managers;
    auto & representation_hypers = Fix::ParentFull::representation_hypers;

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
      for (auto & rep_hypers : representation_hypers[i_manager]) {
        auto & manager = managers[i_manager];
        auto & manager_half = managers_half[i_manager];
        Representation_t representation{rep_hypers};
        representation.compute(manager);
        representation.compute(manager_half);

        auto && rep_vectors{
            *manager->template get_property<Prop_t>(representation.get_name())};
        auto && rep_vectors_half{
            *manager_half->template get_property<PropHalf_t>(
                representation.get_name())};

        auto && rep_vector_gradients{
            *manager->template get_property<PropGrad_t>(
                representation.get_gradient_name())};
        auto && rep_vector_gradients_half{
            *manager_half->template get_property<PropGradHalf_t>(
                representation.get_gradient_name())};

        size_t center_count{0};
        for (auto center : manager) {
          /**
           * To be able to compare the computed values of the full neighbour
           * minimal neighbour list with the values in the full list, a second
           * iterator over the manager with the minimal neighbour list is
           * needed. `center_half` refers to the same `center` atom structure
           * wise, but its values have been computed making using of a minimal
           * neighbourlist.
           */
          auto it_half = manager_half->get_iterator_at(center_count, 0);
          auto center_half = *(it_half);
          // compare the representation coefficients
          auto diff_rep_m{math::relative_error(
              rep_vectors.get_dense_row(center),
              rep_vectors_half.get_dense_row(center_half), delta, epsilon)};
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

          auto ii_pair = center.get_atom_ii();
          auto ii_half_pair = center_half.get_atom_ii();

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

            auto half_neigh_it = center_half.pairs().begin();
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
