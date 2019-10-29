/**
 * @file   test_calculator.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef TESTS_TEST_CALCULATOR_HH_
#define TESTS_TEST_CALCULATOR_HH_

#include "test_adaptor.hh"
#include "test_math.hh"
#include "test_structure.hh"

#include "rascal/atomic_structure.hh"
#include "rascal/json_io.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_covariants.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/cluster_ref_key.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils.hh"

#include <memory>
#include <tuple>

namespace rascal {

  struct TestData {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    TestData() = default;

    void get_ref(const std::string & ref_filename) {
      this->ref_data =
          json::from_ubjson(internal::read_binary_file(ref_filename));
      auto filenames =
          this->ref_data.at("filenames").get<std::vector<std::string>>();
      auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

      for (auto && filename : filenames) {
        for (auto && cutoff : cutoffs) {
          json parameters;
          json structure{{"filename", filename}};
          json adaptors;
          json ad1{{"name", "AdaptorNeighbourList"},
                   {"initialization_arguments", {{"cutoff", cutoff}}}};
          json ad1b{{"name", "AdaptorCenterContribution"},
                    {"initialization_arguments", {}}};
          json ad2{{"name", "AdaptorStrict"},
                   {"initialization_arguments", {{"cutoff", cutoff}}}};
          adaptors.emplace_back(ad1);
          adaptors.emplace_back(ad1b);
          adaptors.emplace_back(ad2);

          parameters["structure"] = structure;
          parameters["adaptors"] = adaptors;

          this->factory_args.emplace_back(parameters);
        }
      }
    }

    ~TestData() = default;

    json ref_data{};
    json factory_args{};
  };

  template <typename MultipleStructureFixture>
  struct MultipleStructureSphericalInvariants : MultipleStructureFixture {
    using Parent = MultipleStructureFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalInvariants;

    MultipleStructureSphericalInvariants() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~MultipleStructureSphericalInvariants() = default;

    std::vector<json> representation_hypers{};

    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};

    std::vector<json> rep_hypers{{{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 3},
                                  {"max_angular", 3},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 4},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 3},
                                  {"max_angular", 1},
                                  {"soap_type", "BiSpectrum"},
                                  {"inversion_symmetry", true},
                                  {"normalize", true}},
                                 {{"max_radial", 3},
                                  {"max_angular", 1},
                                  {"soap_type", "BiSpectrum"},
                                  {"inversion_symmetry", false},
                                  {"normalize", true}}};
  };

  template <typename MultipleStructureFixture>
  struct MultipleStructureSphericalCovariants : MultipleStructureFixture {
    using Parent = MultipleStructureFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalCovariants;

    MultipleStructureSphericalCovariants() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~MultipleStructureSphericalCovariants() = default;

    std::vector<json> representation_hypers{};

    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "AA"}}}},
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{{{"max_radial", 1},
                                  {"max_angular", 2},
                                  {"soap_type", "LambdaSpectrum"},
                                  {"lam", 2},
                                  {"inversion_symmetry", true},
                                  {"normalize", true}},
                                 {{"max_radial", 2},
                                  {"max_angular", 2},
                                  {"soap_type", "LambdaSpectrum"},
                                  {"lam", 2},
                                  {"inversion_symmetry", false},
                                  {"normalize", true}}};
  };

  struct SphericalInvariantsTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalInvariants;
    SphericalInvariantsTestData() : Parent{} {
      this->get_ref(this->ref_filename);
    }
    ~SphericalInvariantsTestData() = default;
    bool verbose{false};
    std::string ref_filename{
        "reference_data/tests_only/spherical_invariants_reference.ubjson"};
  };

  struct SphericalCovariantsTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalCovariants;
    SphericalCovariantsTestData() : Parent{} {
      this->get_ref(this->ref_filename);
    }
    ~SphericalCovariantsTestData() = default;
    bool verbose{false};
    std::string ref_filename{
        "reference_data/tests_only/spherical_covariants_reference.ubjson"};
  };

  template <class MultipleStructureFixture>
  struct MultipleStructureSphericalExpansion : MultipleStructureFixture {
    using Parent = MultipleStructureFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalExpansion;

    MultipleStructureSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~MultipleStructureSphericalExpansion() = default;

    std::vector<json> representation_hypers{};

    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}},
        {{"type", "RadialScaling"},
         {"cutoff", {{"value", 4.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}},
         {"rate", {{"value", .0}, {"unit", "AA"}}},
         {"exponent", {{"value", 4}, {"unit", ""}}},
         {"scale", {{"value", 2.5}, {"unit", "AA"}}}},
        {{"type", "RadialScaling"},
         {"cutoff", {{"value", 4.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}},
         {"rate", {{"value", 1.}, {"unit", "AA"}}},
         {"exponent", {{"value", 3}, {"unit", ""}}},
         {"scale", {{"value", 2.}, {"unit", "AA"}}}}};

    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}},
                                                 {{"type", "DVR"}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.5}, {"unit", "AA"}}}}};

    std::vector<json> rep_hypers{{{"max_radial", 6}, {"max_angular", 4}}};
  };

  /**
   * Simplified version of MultipleStructureManagerNLStrictFixture
   *  that uses only one structure, cutoff, and adaptor set
   *
   *  Useful if we just need a StructureManager to test relatively isolated
   *  functionality on a single structure, but using the rest of the testing
   *  machinery
   */
  struct SimpleStructureManagerNLCCStrictFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;

    SimpleStructureManagerNLCCStrictFixture() {
      json parameters;
      json structure{{"filename", filename}};
      json adaptors;
      json ad1{{"name", "AdaptorNeighbourList"},
               {"initialization_arguments",
                {{"cutoff", cutoff}, {"skin", cutoff_skin}}}};
      json ad1b{{"name", "AdaptorCenterContribution"},
                {"initialization_arguments", {}}};
      json ad2{{"name", "AdaptorStrict"},
               {"initialization_arguments", {{"cutoff", cutoff}}}};
      adaptors.emplace_back(ad1);
      adaptors.push_back(ad1b);
      adaptors.emplace_back(ad2);

      parameters["structure"] = structure;
      parameters["adaptors"] = adaptors;

      this->factory_args.emplace_back(parameters);
    }

    ~SimpleStructureManagerNLCCStrictFixture() = default;

    const std::string filename{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
    const double cutoff{3.};
    const double cutoff_skin{0.};

    json factory_args{};
  };

  struct MultipleHypersSphericalExpansion
      : SimpleStructureManagerNLCCStrictFixture {
    using Parent = SimpleStructureManagerNLCCStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalExpansion;

    MultipleHypersSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~MultipleHypersSphericalExpansion() = default;

    std::vector<json> representation_hypers{};
    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "AA"}}}},
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{
        {{"type", "GTO"}, {"optimization", {{"type", "None"}}}},
        {{"type", "DVR"}, {"optimization", {{"type", "None"}}}},
        {{"type", "GTO"},
         {"optimization",
          {{"type", "Spline"},
           {"accuracy", 1e-12},
           {"range", {{"begin", 0.}, {"end", 3.}}}}}},
        {{"type", "DVR"},
         {"optimization",
          {{"type", "Spline"},
           {"accuracy", 1e-5},
           {"range", {{"begin", 0.000001}, {"end", 3.}}}}}}};
    std::vector<json> rep_hypers{
        {{"max_radial", 4}, {"max_angular", 2}, {"compute_gradients", true}},
        {{"max_radial", 6}, {"max_angular", 4}, {"compute_gradients", true}}};
  };

  /** Contains some simple periodic structures for testing complicated things
   *  like gradients
   */
  struct SimplePeriodicNLCCStrictFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>;
    using Structure_t = AtomicStructure<3>;

    SimplePeriodicNLCCStrictFixture() {
      for (auto && filename : filenames) {
        json parameters;
        json structure{{"filename", filename}};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                 {"initialization_arguments",
                  {{"cutoff", cutoff}, {"skin", cutoff_skin}}}};
        json ad1b{{"name", "AdaptorCenterContribution"},
                  {"initialization_arguments", {}}};
        json ad2{{"name", "AdaptorStrict"},
                 {"initialization_arguments", {{"cutoff", cutoff}}}};
        adaptors.emplace_back(ad1);
        adaptors.push_back(ad1b);
        adaptors.emplace_back(ad2);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);
      }
    }

    ~SimplePeriodicNLCCStrictFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/simple_cubic_8.json",
        "reference_data/diamond_2atom_distorted.json",
        "reference_data/diamond_cubic_distorted.json",
        "reference_data/SiCGe_wurtzite_like.json",
        "reference_data/SiC_moissanite_supercell.json",
        "reference_data/small_molecule.json",
        "reference_data/methane.json"};
    // Simpler structures for debugging:
    // "reference_data/diamond_2atom.json",
    // "reference_data/SiC_moissanite.json",
    const double cutoff{2.5};
    const double cutoff_skin{0.};

    json factory_args{};
    std::vector<Structure_t> structures{};
  };

  struct SingleHypersSphericalExpansion : SimplePeriodicNLCCStrictFixture {
    using Parent = SimplePeriodicNLCCStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalExpansion;

    SingleHypersSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~SingleHypersSphericalExpansion() = default;

    std::vector<json> representation_hypers{};
    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 2.5}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{
        {{"max_radial", 2}, {"max_angular", 2}, {"compute_gradients", true}},
        {{"max_radial", 3}, {"max_angular", 0}, {"compute_gradients", true}}};
  };

  struct SingleHypersSphericalInvariants : SimplePeriodicNLCCStrictFixture {
    using Parent = SimplePeriodicNLCCStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalInvariants;

    SingleHypersSphericalInvariants() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~SingleHypersSphericalInvariants() = default;

    std::vector<json> representation_hypers{};
    std::vector<json> fc_hypers{
        {{"type", "ShiftedCosine"},
         {"cutoff", {{"value", 2.5}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

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

  struct SphericalExpansionTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalExpansion;

    SphericalExpansionTestData() : Parent{} {
      this->get_ref(this->ref_filename);
    }
    ~SphericalExpansionTestData() = default;
    bool verbose{false};
    std::string ref_filename{
        "reference_data/tests_only/spherical_expansion_reference.ubjson"};
  };

  /**
   * Calculator specialized to testing the derivative of the RadialIntegral
   * in the definition of the SphericalExpansion representation.
   */
  template <class RadialIntegral, class ClusterRef>
  struct SphericalExpansionRadialDerivative {
    SphericalExpansionRadialDerivative(std::shared_ptr<RadialIntegral> ri,
                                       ClusterRef & pair_in)
        : radial_integral{ri}, pair{pair_in}, max_radial{ri->max_radial},
          max_angular{ri->max_angular} {}

    ~SphericalExpansionRadialDerivative() = default;

    Eigen::Array<double, 1, Eigen::Dynamic>
    f(const Eigen::Matrix<double, 1, 1> & input_v) {
      Eigen::ArrayXXd result(this->max_radial, this->max_angular + 1);
      result = this->radial_integral->template compute_neighbour_contribution(
          input_v(0), this->pair);
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    Eigen::Array<double, 1, Eigen::Dynamic>
    grad_f(const Eigen::Matrix<double, 1, 1> & input_v) {
      Eigen::ArrayXXd result(this->max_radial, this->max_angular + 1);
      result = this->radial_integral->template compute_neighbour_derivative(
          input_v(0), this->pair);
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    std::shared_ptr<RadialIntegral> radial_integral;
    ClusterRef & pair;
    size_t max_radial{6};
    size_t max_angular{4};
  };

  template <class BaseFixture, internal::RadialBasisType RadialType,
            internal::AtomicSmearingType SmearingType,
            internal::OptimizationType OptType>
  struct RadialIntegralHandlerFixture : MultipleStructureFixture<BaseFixture> {
    using Parent = MultipleStructureFixture<BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using RadialIntegral_t =
        internal::RadialContributionHandler<RadialType, SmearingType, OptType>;

    RadialIntegralHandlerFixture() : Parent{} {
      // filter out the hypers that don't correspond to the current RadialType,
      // SmearingType or OptType
      std::vector<json> hypers_temp;

      for (const auto & hyper : this->representation_hypers) {
        // This block is to ignore hypers which do not agree with the type of
        // the fixture. This way we do not have to create a fixture for each
        // type while not using the wrong templated RadialIntegralHandler
        // constructor
        auto radial_contribution_hypers =
            hyper.at("radial_contribution").template get<json>();
        auto radial_contribution_type_name =
            radial_contribution_hypers.at("type").template get<std::string>();
        auto density_hypers = hyper.at("gaussian_density").template get<json>();
        auto smearing_type_name =
            density_hypers.at("type").template get<std::string>();
        auto optimization_hypers =
            radial_contribution_hypers.at("optimization").template get<json>();
        auto optimization_type_name =
            optimization_hypers.at("type").template get<std::string>();

        internal::RadialBasisType radial_contribution_type{};
        internal::AtomicSmearingType smearing_type{};
        internal::OptimizationType optimization_type{};
        if (radial_contribution_type_name == "GTO") {
          radial_contribution_type = internal::RadialBasisType::GTO;
        } else if (radial_contribution_type_name == "DVR") {
          radial_contribution_type = internal::RadialBasisType::DVR;
        } else {
          throw std::runtime_error(
              "Wrong radial basis type for RadialIntegralHandler tests");
        }

        if (smearing_type_name == "Constant") {
          smearing_type = internal::AtomicSmearingType::Constant;
        } else {
          throw std::runtime_error(
              "Wrong smearing type for RadialIntegralHandler tests");
        }

        if (optimization_type_name == "None") {
          optimization_type = internal::OptimizationType::None;
        } else if (optimization_type_name == "Spline") {
          optimization_type = internal::OptimizationType::Interpolator;
        } else {
          throw std::runtime_error(
              "Wrong optimization type for RadialIntegralHandler tests");
        }
        auto hypers_radial_contribution_handler_type{
            internal::combine_to_radial_contribution_type(
                radial_contribution_type, smearing_type, optimization_type)};

        auto radial_contribution_handler_type{
            internal::combine_to_radial_contribution_type(
                RadialType, SmearingType, OptType)};

        if (hypers_radial_contribution_handler_type ==
            radial_contribution_handler_type) {
          hypers_temp.push_back(hyper);
        }
      }
      this->representation_hypers.clear();
      this->representation_hypers = std::move(hypers_temp);
    }
    ~RadialIntegralHandlerFixture() = default;
  };

  /**
   * Gradient provider specialized to testing the gradient of a Calculator
   *
   * The gradient is tested center-by-center, by iterating over each center and
   * doing finite displacements on its position.  This iteration should normally
   * be done by the RepresentationCalculatorGradientFixture class.
   *
   * In the case of periodic structures, the gradient is accumulated only onto
   * _real_ atoms, but the motion of all _images_ of the "moving" atom (the one
   * with respect to which the gradient is being taken) is taken into account.
   *
   * Initialize with a Calculator, a StructureManager, and an
   * AtomicStructure representing the original structure (before modifying with
   * finite-difference displacements).  The gradient of the representation with
   * respect to the center position can then be tested, as usual, with
   * test_gradients() (defined in test_math.hh).
   */
  template <typename Calculator, class StructureManager>
  class RepresentationCalculatorGradientProvider {
   public:
    using Structure_t = AtomicStructure<3>;
    using Key_t = typename Calculator::Key_t;
    static const size_t n_arguments = 3;

    using PairRef_t =
        typename Calculator::template ClusterRef_t<StructureManager, 2>;

    using PairRefKey_t = typename PairRef_t::ThisParentClass;

    // type of the data structure holding the representation and its gradients
    using Prop_t = typename Calculator::template Property_t<StructureManager>;
    using PropGrad_t =
        typename Calculator::template PropertyGradient_t<StructureManager>;

    template <typename T, class V>
    friend class RepresentationCalculatorGradientFixture;

    RepresentationCalculatorGradientProvider(
        Calculator & representation,
        std::shared_ptr<StructureManager> structure_manager,
        Structure_t atomic_structure)
        : representation{representation}, structure_manager{structure_manager},
          atomic_structure{atomic_structure}, center_it{
                                                  structure_manager->begin()} {
      for (auto center : this->structure_manager) {
        this->n_neighbors.push_back(center.size());
      }
    }

    ~RepresentationCalculatorGradientProvider() = default;

    Eigen::Array<double, 1, Eigen::Dynamic>
    f(const Eigen::Ref<const Eigen::Vector3d> & center_position) {
      auto center = *center_it;
      Structure_t modified_structure{this->atomic_structure};
      modified_structure.positions.col(center.get_index()) = center_position;
      modified_structure.wrap();
      this->structure_manager->update(modified_structure);
      int i_center{0};
      for (auto center : this->structure_manager) {
        if (this->n_neighbors[i_center] != center.size()) {
          throw std::runtime_error(
              R"(The number of neighbors has changed when making finite
              displacements. This happens because a neighbor is almost at the
              cutoff boundary so please change the structure or the cutoff to
              avoid this.)");
        }
        ++i_center;
      }

      this->representation.compute(this->structure_manager);

      auto && data_sparse{*structure_manager->template get_property_ptr<Prop_t>(
          representation.get_name())};
      auto && gradients_sparse{
          *structure_manager->template get_property_ptr<PropGrad_t>(
              representation.get_gradient_name())};
      auto ii_pair = center.get_atom_ii();
      auto & data_center{data_sparse[ii_pair]};
      auto keys_center = gradients_sparse.get_keys(ii_pair);
      Key_t center_key{center.get_atom_type()};
      size_t n_entries_per_key{static_cast<size_t>(data_sparse.get_nb_comp())};
      size_t n_entries_center{n_entries_per_key * keys_center.size()};
      size_t n_entries_neighbours{0};
      // Count all the keys in the sparse gradient structure where the gradient
      // is nonzero (i.e. where the key has an entry in the structure)
      for (auto neigh : center) {
        if (this->structure_manager->is_ghost_atom(neigh)) {
          // Don't compute gradient contributions onto ghost atoms
          continue;
        }
        auto swapped_ref{std::move(swap_pair_ref(neigh).front())};
        n_entries_neighbours +=
            (gradients_sparse[swapped_ref].get_keys().size() *
             n_entries_per_key);
      }
      // Packed array containing: The center coefficients (all species) and
      // the neighbour coefficients (only same species as center)
      Eigen::ArrayXd data_pairs(n_entries_center + n_entries_neighbours);

      size_t result_idx{0};
      for (auto & key : keys_center) {
        Eigen::Map<Eigen::RowVectorXd> data_flat(data_center[key].data(),
                                                 n_entries_per_key);
        data_pairs.segment(result_idx, n_entries_per_key) = data_flat;
        result_idx += n_entries_per_key;
      }
      for (auto neigh : center) {
        if (this->structure_manager->is_ghost_atom(neigh)) {
          // Don't compute gradient contributions onto ghost atoms
          continue;
        }
        auto & data_neigh{data_sparse[neigh]};
        // The neighbour gradient (i =/= j) only contributes to certain species
        // channels (keys), in the case of SOAP and SphExpn those keys
        // containing the species of the center (the atom wrt the derivative is
        // being taken)
        // The nonzero gradient keys are already indicated in the sparse
        // gradient structure
        auto swapped_ref{std::move(swap_pair_ref(neigh).front())};
        auto keys_neigh{gradients_sparse[swapped_ref].get_keys()};
        for (auto & key : keys_neigh) {
          Eigen::Map<Eigen::ArrayXd> data_flat(data_neigh[key].data(),
                                               n_entries_per_key);
          data_pairs.segment(result_idx, n_entries_per_key) = data_flat;
          result_idx += n_entries_per_key;
        }
      }

      // Reset the atomic structure for the next iteration
      this->structure_manager->update(this->atomic_structure);
      return data_pairs.transpose();
    }

    Eigen::Array<double, 3, Eigen::Dynamic>
    grad_f(const Eigen::Ref<const Eigen::Vector3d> & /*center_position*/) {
      using Matrix3Xd_RowMaj_t =
          Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
      // Assume f() was already called and updated the position
      // center_it->position() = center_position;
      // representation.compute();
      auto center = *center_it;

      auto && data_sparse{*structure_manager->template get_property_ptr<Prop_t>(
          representation.get_name())};
      auto && gradients_sparse{
          *structure_manager->template get_property_ptr<PropGrad_t>(
              representation.get_gradient_name())};
      auto ii_pair = center.get_atom_ii();
      auto & gradients_center{gradients_sparse[ii_pair]};
      auto keys_center = gradients_center.get_keys();
      size_t n_entries_per_key{static_cast<size_t>(data_sparse.get_nb_comp())};
      size_t n_entries_center{n_entries_per_key * keys_center.size()};
      size_t n_entries_neighbours{0};
      for (auto neigh : center) {
        if (this->structure_manager->is_ghost_atom(neigh)) {
          // Don't compute gradient contributions onto ghost atoms
          continue;
        }
        auto swapped_ref{std::move(swap_pair_ref(neigh).front())};
        n_entries_neighbours +=
            (gradients_sparse[swapped_ref].get_keys().size() *
             n_entries_per_key);
      }
      Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>
          grad_coeffs_pairs(3, n_entries_center + n_entries_neighbours);
      grad_coeffs_pairs.setZero();

      // Use the exact same iteration pattern as in f()  to guarantee that the
      // gradients appear in the same place as their corresponding data
      size_t result_idx{0};
      for (auto & key : keys_center) {
        // Here the 'flattening' retains the 3 Cartesian dimensions as rows,
        // since they vary the slowest within each key
        Eigen::Map<Matrix3Xd_RowMaj_t> grad_coeffs_flat(
            gradients_center[key].data(), 3, n_entries_per_key);
        grad_coeffs_pairs.block(0, result_idx, 3, n_entries_per_key) =
            grad_coeffs_flat;
        result_idx += n_entries_per_key;
      }
      for (auto neigh : center) {
        if (this->structure_manager->is_ghost_atom(neigh)) {
          // Don't compute gradient contributions onto ghost atoms
          continue;
        }
        // We need grad_i c^{ji} -- using just 'neigh' would give us
        // grad_j c^{ij}, hence the swap
        auto neigh_swap_images{swap_pair_ref(neigh)};
        auto & gradients_neigh_first{
            gradients_sparse[neigh_swap_images.front()]};
        // The set of species keys should be the same for all images of i
        auto keys_neigh{gradients_neigh_first.get_keys()};
        for (auto & key : keys_neigh) {
          // For each key, accumulate gradients over periodic images of the atom
          // that moves in the finite-difference step
          for (auto & neigh_swap : neigh_swap_images) {
            auto & gradients_neigh{gradients_sparse[neigh_swap]};
            Eigen::Map<Matrix3Xd_RowMaj_t> grad_coeffs_flat(
                gradients_neigh[key].data(), 3, n_entries_per_key);
            grad_coeffs_pairs.block(0, result_idx, 3, n_entries_per_key) +=
                grad_coeffs_flat;
          }
          result_idx += n_entries_per_key;
        }
      }
      return grad_coeffs_pairs;
    }

   private:
    Calculator & representation;
    std::shared_ptr<StructureManager> structure_manager;
    Structure_t atomic_structure;
    typename StructureManager::iterator center_it;
    //! count the number of neighbours of each centers
    std::vector<size_t> n_neighbors{};

    void advance_center() { ++this->center_it; }

    /**
     * Swap a ClusterRef<order=2> (i, j) so it refers to (j, i) instead
     *
     * @return std::vector of ClusterRefKeys or order 2 (pair keys) of all pairs
     *         (j, i') where i' is either i or any of its periodic images within
     *         the cutoff of j. The atom j, on the other hand, must be a real
     *         atom (not a ghost or periodic image).
     */
    std::vector<PairRefKey_t> swap_pair_ref(const PairRef_t & pair_ref) {
      auto center_manager{extract_underlying_manager<0>(structure_manager)};
      auto atomic_structure{center_manager->get_atomic_structure()};
      // Get the atom index to the corresponding atom tag
      size_t access_index{structure_manager->get_atom_index(pair_ref.back())};
      auto new_center_it{structure_manager->get_iterator_at(access_index)};
      // Return cluster ref at which the iterator is currently pointing
      auto && new_center{*new_center_it};
      size_t i_index{structure_manager->get_atom_index(pair_ref.front())};

      // Find all (j, i') pairs
      std::vector<PairRefKey_t> new_pairs;
      for (auto new_pair : new_center) {
        size_t i_trial_index{
            structure_manager->get_atom_index(new_pair.back())};
        // Is this the i (old center) atom or any of its images?
        if (i_trial_index == i_index) {
          new_pairs.emplace_back(std::move(new_pair));
        }
      }
      if (new_pairs.size() == 0) {
        std::stringstream err_str{};
        err_str << "Didn't find any pairs for pair (i=" << pair_ref.front()
                << ", j=" << pair_ref.back()
                << "); access index for j = " << access_index;
        throw std::range_error(err_str.str());
      }
      return new_pairs;
    }
  };

  /**
   * Test fixture holding the gradient calculator and structure manager
   *
   * Holds data (function values, gradient directions, verbosity) and iterates
   * through the list of centers
   */
  template <typename Calculator, class StructureManager>
  class RepresentationCalculatorGradientFixture : public GradientTestFixture {
   public:
    using StdVector2Dim_t = std::vector<std::vector<double>>;
    using Provider_t =
        RepresentationCalculatorGradientProvider<Calculator, StructureManager>;

    static const size_t n_arguments = 3;
    /**
     * Increased error tolerance because some representations have quite large
     * finite-difference truncation errors (and possibly numerical issues for
     * very small displacements)
     */
    double fd_error_tol{1E-4};

    /**
     * Initialize a gradient test fixture
     *
     * @param filename JSON file holding gradient test parameters, format
     *                 documented in GradientTestFixture
     *
     * @param structure StructureManager on which to test
     *
     * @param calc RepresentationCalculator whose gradient is being tested
     */
    RepresentationCalculatorGradientFixture(
        std::string filename, std::shared_ptr<StructureManager> structure,
        Provider_t & calc)
        : structure{structure}, center_it{structure->begin()}, provider{calc} {
      json input_data = json_io::load(filename);

      this->function_inputs = this->get_function_inputs();
      this->displacement_directions =
          this->get_displacement_directions(input_data, this->n_arguments);
      this->verbosity = get_verbosity(input_data);
      if (input_data.find("fd_error_tol") != input_data.end()) {
        this->fd_error_tol = input_data["fd_error_tol"].get<double>();
      }
    }

    ~RepresentationCalculatorGradientFixture() = default;

    const Provider_t & get_provider() { return provider; }

    /**
     * Go to the next center in the structure
     *
     * Not (yet) implemented as iterator because that overcomplicates things
     */
    void advance_center() {
      ++this->center_it;
      this->provider.advance_center();
      if (this->has_next()) {
        this->function_inputs = get_function_inputs();
      }
    }

    bool has_next() { return (this->center_it != structure->end()); }

   private:
    StdVector2Dim_t get_function_inputs() {
      StdVector2Dim_t inputs_new{};
      auto center_pos = (*center_it).get_position();
      inputs_new.emplace_back(center_pos.data(),
                              center_pos.data() + center_pos.size());
      return inputs_new;
    }

    std::shared_ptr<StructureManager> structure;
    typename StructureManager::iterator center_it;
    Provider_t & provider;
  };

  template <class MultipleStructureFixture>
  struct MultipleStructureSortedCoulomb : MultipleStructureFixture {
    using Parent = MultipleStructureFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSortedCoulomb;
    MultipleStructureSortedCoulomb() : Parent{} {};
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<json> representation_hypers{
        {{"central_cutoff", 3.},
         {"central_decay", 0.5},
         {"interaction_cutoff", 10.},
         {"interaction_decay", 0.5},
         {"size", 120},
         {"sorting_algorithm", "distance"}},
        {{"central_cutoff", 3.},
         {"central_decay", 0.5},
         {"interaction_cutoff", 10.},
         {"interaction_decay", 0.5},
         {"size", 120},
         {"sorting_algorithm", "row_norm"}}};
  };

  struct SortedCoulombTestData {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    using Representation_t = CalculatorSortedCoulomb;
    SortedCoulombTestData() { this->get_ref(this->ref_filename); }
    ~SortedCoulombTestData() = default;

    void get_ref(const std::string & ref_filename) {
      this->ref_data =
          json::from_ubjson(internal::read_binary_file(ref_filename));
      auto filenames =
          this->ref_data.at("filenames").get<std::vector<std::string>>();
      auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

      for (auto && filename : filenames) {
        for (auto && cutoff : cutoffs) {
          // std::cout << filename << " " << cutoff << std::endl;
          json parameters;
          json structure{{"filename", filename}};
          json adaptors;
          json ad1{{"name", "AdaptorNeighbourList"},
                   {"initialization_arguments", {{"cutoff", cutoff}}}};
          json ad2{{"name", "AdaptorStrict"},
                   {"initialization_arguments", {{"cutoff", cutoff}}}};
          adaptors.emplace_back(ad1);
          adaptors.emplace_back(ad2);

          parameters["structure"] = structure;
          parameters["adaptors"] = adaptors;

          this->factory_args.emplace_back(parameters);
        }
      }
    }

    json ref_data{};
    json factory_args{};

    std::string ref_filename{
        "reference_data/tests_only/sorted_coulomb_reference.ubjson"};
    bool verbose{false};
  };

  template <class BaseFixture>
  struct CalculatorFixture : MultipleStructureFixture<BaseFixture> {
    using Parent = MultipleStructureFixture<BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = typename BaseFixture::Representation_t;
    using Property_t =
        typename Representation_t::template Property_t<Manager_t>;

    CalculatorFixture() : Parent{} {}
    ~CalculatorFixture() = default;

    std::vector<Representation_t> representations{};
  };

  /* ---------------------------------------------------------------------- */

}  // namespace rascal

#endif  // TESTS_TEST_CALCULATOR_HH_
