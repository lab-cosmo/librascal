/**
 * file   test_calculator.hh
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

#include "tests.hh"
#include "test_adaptor.hh"
#include "test_math.hh"
#include "test_structure.hh"
#include "atomic_structure.hh"
#include "structure_managers/structure_manager_collection.hh"
#include "representations/calculator_base.hh"
#include "representations/calculator_sorted_coulomb.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"
#include "representations/calculator_spherical_covariants.hh"

#include "json_io.hh"
#include "rascal_utility.hh"

#include <tuple>
#include <memory>

namespace rascal {

  struct TestData {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    TestData() = default;

    void get_ref(const std::string & ref_filename) {
      std::vector<std::uint8_t> ref_data_ubjson;
      internal::read_binary_file(ref_filename, ref_data_ubjson);
      this->ref_data = json::from_ubjson(ref_data_ubjson);
      auto filenames =
          this->ref_data.at("filenames").get<std::vector<std::string>>();
      auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

      for (auto && filename : filenames) {
        for (auto && cutoff : cutoffs) {
          // std::cout << filename << " " << cutoff << std::endl;
          json parameters;
          json structure{{"filename", filename}};
          json adaptors;
          json ad1{
              {"name", "AdaptorNeighbourList"},
              {"initialization_arguments",
               {{"cutoff", cutoff},
                {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
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

    ~TestData() = default;

    const bool consider_ghost_neighbours{false};
    json ref_data{};
    json factory_args{};
  };

  struct MultipleStructureSphericalInvariants
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
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
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "AA"}}}},
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};

    std::vector<json> rep_hypers{{{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 10},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 3},
                                  {"max_angular", 3},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 8},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 4},
                                  {"max_angular", 1},
                                  {"soap_type", "BiSpectrum"},
                                  {"inversion_symmetry", true},
                                  {"normalize", true}},
                                 {{"max_radial", 4},
                                  {"max_angular", 1},
                                  {"soap_type", "BiSpectrum"},
                                  {"inversion_symmetry", false},
                                  {"normalize", true}}};
  };

  struct MultipleStructureSphericalCovariants
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
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
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "Cosine"},
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
                                  {"max_angular", 3},
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
        "reference_data/spherical_invariants_reference.ubjson"};
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
        "reference_data/spherical_covariants_reference.ubjson"};
  };

  struct MultipleStructureSphericalExpansion
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
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
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.5}, {"unit", "AA"}}}}};

    std::vector<json> rep_hypers{{{"max_radial", 10}, {"max_angular", 8}}};
  };

  /** Simplified version of MultipleStructureManagerNLStrictFixture
   *  that uses only one structure, cutoff, and adaptor set
   *
   *  Useful if we just need a StructureManager to test relatively isolated
   *  functionality on a single structure, but using the rest of the testing
   *  machinery
   */
  struct SimpleStructureManagerNLStrictFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;

    SimpleStructureManagerNLStrictFixture() {
      json parameters;
      json structure{{"filename", filename}};
      json adaptors;
      json ad1{{"name", "AdaptorNeighbourList"},
               {"initialization_arguments",
                {{"cutoff", cutoff},
                 {"skin", cutoff_skin},
                 {"consider_ghost_neighbours", false}}}};
      json ad2{{"name", "AdaptorStrict"},
               {"initialization_arguments", {{"cutoff", cutoff}}}};
      adaptors.emplace_back(ad1);
      adaptors.emplace_back(ad2);

      parameters["structure"] = structure;
      parameters["adaptors"] = adaptors;

      this->factory_args.emplace_back(parameters);
    }

    ~SimpleStructureManagerNLStrictFixture() = default;

    const std::string filename{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
    const double cutoff{3.};
    const double cutoff_skin{0.5};

    json factory_args{};
  };

  struct MultipleHypersSphericalExpansion
      : SimpleStructureManagerNLStrictFixture {
    using Parent = SimpleStructureManagerNLStrictFixture;
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
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}},
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "AA"}}}},
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{
        {{"max_radial", 4}, {"max_angular", 2}, {"compute_gradients", true}},
        {{"max_radial", 6}, {"max_angular", 4}, {"compute_gradients", true}}};
  };

  /** Contains some simple periodic structures for testing complicated things
   *  like gradients
   */
  struct SimplePeriodicNLStrictFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    using Structure_t = AtomicStructure<3>;

    SimplePeriodicNLStrictFixture() {
      for (auto && filename : filenames) {
        json parameters;
        json structure{{"filename", filename}};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                 {"initialization_arguments",
                  {{"cutoff", cutoff},
                   {"skin", cutoff_skin},
                   {"consider_ghost_neighbours", false}}}};
        json ad2{{"name", "AdaptorStrict"},
                 {"initialization_arguments", {{"cutoff", cutoff}}}};
        adaptors.emplace_back(ad1);
        adaptors.emplace_back(ad2);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);
      }
    }

    ~SimplePeriodicNLStrictFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/diamond_2atom.json",
        "reference_data/diamond_2atom_distorted.json",
        "reference_data/diamond_cubic_distorted.json",
        "reference_data/SiC_moissanite.json",
        "reference_data/SiCGe_wurtzite_like.json",
        "reference_data/SiC_moissanite_supercell.json"};
    const double cutoff{2.5};
    const double cutoff_skin{0.5};

    json factory_args{};
    std::vector<Structure_t> structures{};
  };

  struct SingleHypersSphericalExpansion : SimplePeriodicNLStrictFixture {
    using Parent = SimplePeriodicNLStrictFixture;
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
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.5}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{
        {{"max_radial", 2}, {"max_angular", 2}, {"compute_gradients", true}},
        {{"max_radial", 4}, {"max_angular", 0}, {"compute_gradients", true}}};
  };

  struct SingleHypersSphericalInvariants : SimplePeriodicNLStrictFixture {
    using Parent = SimplePeriodicNLStrictFixture;
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
        {{"type", "Cosine"},
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
                                 {{"max_radial", 4},
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
        "reference_data/spherical_expansion_reference.ubjson"};
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
      result = this->radial_integral->template compute_neighbour_contribution<
          internal::AtomicSmearingType::Constant>(input_v(0), this->pair);
      // result.matrix().transpose() *=
      //     this->radial_integral->radial_norm_factors.asDiagonal();
      // result.matrix().transpose() *=
      //     this->radial_integral->radial_ortho_matrix;
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    Eigen::Array<double, 1, Eigen::Dynamic>
    grad_f(const Eigen::Matrix<double, 1, 1> & input_v) {
      Eigen::ArrayXXd result(this->max_radial, this->max_angular + 1);
      result = this->radial_integral->template compute_neighbour_derivative<
          internal::AtomicSmearingType::Constant>(input_v(0), this->pair);
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    std::shared_ptr<RadialIntegral> radial_integral;
    ClusterRef & pair;
    size_t max_radial{6};
    size_t max_angular{4};
  };

  /**
   * Calculator specialized to testing the gradient of a RepresentationManager
   *
   * The gradient is tested center-by-center, by iterating over each center and
   * doing finite displacements on its position.  This iteration should normally
   * be done by the RepresentationManagerGradientFixture class.
   *
   * Initialize with a RepresentationManager, a StructureManager, and an
   * AtomicStructure representing the original structure (before modifying with
   * finite-difference displacments).  The gradient of the representation with
   * respect to the center position can then be tested, as usual, with
   * test_gradients() (defined in test_math.hh).
   */
  template <typename RepManager, class StructureManager>
  class RepresentationManagerGradientCalculator {
   public:
    using Structure_t = AtomicStructure<3>;
    using Key_t = typename RepManager::Key_t;
    static const size_t n_arguments = 3;

    using PairRef_t =
        typename RepManager::template ClusterRef_t<StructureManager, 2>;

    // type of the data structure holding the representation and its gradients
    using Prop_t = typename RepManager::template Property_t<StructureManager>;
    using PropGrad_t =
        typename RepManager::template PropertyGradient_t<StructureManager>;

    template <typename T, class V>
    friend class RepresentationManagerGradientFixture;

    RepresentationManagerGradientCalculator(
        RepManager & representation,
        std::shared_ptr<StructureManager> structure_manager,
        Structure_t atomic_structure)
        : representation{representation}, structure_manager{structure_manager},
          atomic_structure{atomic_structure}, center_it{
                                                  structure_manager->begin()} {}

    ~RepresentationManagerGradientCalculator() = default;

    Eigen::Array<double, 1, Eigen::Dynamic>
    f(const Eigen::Ref<const Eigen::Vector3d> & center_position) {
      auto center = *center_it;
      Structure_t modified_structure{this->atomic_structure};
      modified_structure.positions.col(center.get_index()) = center_position;
      this->structure_manager->update(modified_structure);
      this->representation.compute(this->structure_manager);

      auto && data_sparse{structure_manager->template get_property_ref<Prop_t>(
          representation.get_name())};
      auto && gradients_sparse{
          structure_manager->template get_property_ref<PropGrad_t>(
              representation.get_gradient_name())};

      auto & data_center{data_sparse[center]};
      auto keys_center = gradients_sparse.get_keys(center);
      Key_t center_key{center.get_atom_type()};
      size_t n_entries_per_key{static_cast<size_t>(data_sparse.get_nb_comp())};
      size_t n_entries_center{n_entries_per_key * keys_center.size()};
      size_t n_entries_neighbours{0};
      // Count all the keys in the sparse gradient structure where the gradient
      // is nonzero (i.e. where the key has an entry in the structure)
      for (auto neigh : center) {
        n_entries_neighbours +=
            (gradients_sparse[swap_pair_ref(neigh)].get_keys().size() *
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
        auto & data_neigh{data_sparse[neigh]};
        // The neighbour gradient (i =/= j) only contributes to certain species
        // channels (keys), in the case of SOAP and SphExpn those keys
        // containing the species of the center (the atom wrt the derivative is
        // being taken)
        // The nonzero gradient keys are already indicated in the sparse
        // gradient structure
        auto keys_neigh{gradients_sparse[swap_pair_ref(neigh)].get_keys()};
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

      auto && data_sparse{structure_manager->template get_property_ref<Prop_t>(
          representation.get_name())};
      auto && gradients_sparse{
          structure_manager->template get_property_ref<PropGrad_t>(
              representation.get_gradient_name())};

      auto & gradients_center{gradients_sparse[center]};
      auto keys_center = gradients_center.get_keys();
      size_t n_entries_per_key{static_cast<size_t>(data_sparse.get_nb_comp())};
      size_t n_entries_center{n_entries_per_key * keys_center.size()};
      size_t n_entries_neighbours{0};
      for (auto neigh : center) {
        n_entries_neighbours +=
            (gradients_sparse[swap_pair_ref(neigh)].get_keys().size() *
             n_entries_per_key);
      }
      Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>
          grad_coeffs_pairs(3, n_entries_center + n_entries_neighbours);

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
        // We need grad_i c^{ji} -- using just 'neigh' would give us
        // grad_j c^{ij}, hence the swap
        auto neigh_swap{swap_pair_ref(neigh)};
        auto & gradients_neigh{gradients_sparse[neigh_swap]};
        auto keys_neigh{gradients_neigh.get_keys()};
        for (auto key : keys_neigh) {
          Eigen::Map<Matrix3Xd_RowMaj_t> grad_coeffs_flat(
              gradients_neigh[key].data(), 3, n_entries_per_key);
          grad_coeffs_pairs.block(0, result_idx, 3, n_entries_per_key) =
              grad_coeffs_flat;
          result_idx += n_entries_per_key;
        }
      }
      return grad_coeffs_pairs;
    }

   private:
    RepManager & representation;
    std::shared_ptr<StructureManager> structure_manager;
    Structure_t atomic_structure;
    typename StructureManager::iterator center_it;

    inline void advance_center() { ++this->center_it; }

    /**
     * Swap a ClusterRef (i, j) so it refers to (j, i) instead
     *
     * @todo wouldn't this be better as a member of StructureManager
     *       (viz. AdaptorNeighbourList<whatever>)?
     */
    PairRef_t swap_pair_ref(const PairRef_t & pair_ref) {
      // Get the atom index to the corresponding atom tag
      size_t access_index = structure_manager->get_atom_index(pair_ref.back());
      auto new_center_it{structure_manager->get_iterator_at(access_index)};
      // Return cluster ref at which the iterator is currently pointing
      auto && new_center{*new_center_it};
      // Iterate until (j,i) is found
      for (auto new_pair : new_center) {
        if (new_pair.back() == pair_ref.front()) {
          return new_pair;
        }
      }
      std::stringstream err_str{};
      err_str << "Didn't find symmetric pair for pair (i=" << pair_ref.front()
              << ", j=" << pair_ref.back() << ").";
      throw std::range_error(err_str.str());
    }
  };

  /**
   * Test fixture holding the gradient calculator and structure manager
   *
   * Holds data (i.e. function values, gradient directions) and iterates through
   * the list of centers
   */
  template <typename RepManager_t, class StructureManager_t>
  class RepresentationManagerGradientFixture : public GradientTestFixture {
   public:
    using StdVector2Dim_t = std::vector<std::vector<double>>;
    using Calculator_t =
        RepresentationManagerGradientCalculator<RepManager_t,
                                                StructureManager_t>;

    static const size_t n_arguments = 3;

    RepresentationManagerGradientFixture(
        std::string filename, std::shared_ptr<StructureManager_t> structure,
        Calculator_t & calc)
        : structure{structure}, center_it{structure->begin()}, calculator{
                                                                   calc} {
      json input_data = json_io::load(filename);

      this->function_inputs = this->get_function_inputs();
      this->displacement_directions =
          this->get_displacement_directions(input_data, this->n_arguments);
      this->verbosity = get_verbosity(input_data);
      if (input_data.find("fd_error_tol") != input_data.end()) {
        this->fd_error_tol = input_data["fd_error_tol"].get<double>();
      }
    }

    ~RepresentationManagerGradientFixture() = default;

    const Calculator_t & get_calculator() { return calculator; }

    /**
     * Go to the next center in the structure
     *
     * Not (yet) implemented as iterator because that overcomplicates things
     */
    inline void advance_center() {
      ++this->center_it;
      this->calculator.advance_center();
      if (this->has_next()) {
        this->function_inputs = get_function_inputs();
      }
    }

    inline bool has_next() { return (this->center_it != structure->end()); }

   private:
    StdVector2Dim_t get_function_inputs() {
      StdVector2Dim_t inputs_new{};
      auto center_pos = (*center_it).get_position();
      inputs_new.emplace_back(center_pos.data(),
                              center_pos.data() + center_pos.size());
      return inputs_new;
    }

    std::shared_ptr<StructureManager_t> structure;
    typename StructureManager_t::iterator center_it;
    Calculator_t & calculator;
  };

  struct MultipleStructureSortedCoulomb
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
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

  struct SortedCoulombTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSortedCoulomb;
    SortedCoulombTestData() : Parent{} { this->get_ref(this->ref_filename); }
    ~SortedCoulombTestData() = default;

    // name of the file containing the reference data. it has been generated
    // with the following python code:
    // script/generate_sorted_coulomb_ref_data.py

    const bool consider_ghost_neighbours{false};
    std::string ref_filename{"reference_data/sorted_coulomb_reference.ubjson"};
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
