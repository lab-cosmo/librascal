/**
 * file   test_representation_manager_base.hh
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

#ifndef TESTS_TEST_REPRESENTATION_MANAGER_HH_
#define TESTS_TEST_REPRESENTATION_MANAGER_HH_

#include "tests.hh"
#include "test_adaptor.hh"
#include "test_math.hh"
#include "test_structure.hh"
#include "atomic_structure.hh"
#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap.hh"
#include "representations/feature_manager_block_sparse.hh"

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

  struct MultipleStructureSOAP : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSOAP() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~MultipleStructureSOAP() = default;

    std::vector<json> hypers{};

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
                                 {{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 1},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}}};
  };

  struct SOAPTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    SOAPTestData() : Parent{} { this->get_ref(this->ref_filename); }
    ~SOAPTestData() = default;
    std::string ref_filename{"reference_data/soap_reference.ubjson"};
  };

  struct MultipleStructureSphericalExpansion
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~MultipleStructureSphericalExpansion() = default;

    std::vector<json> hypers{};

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

    MultipleHypersSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~MultipleHypersSphericalExpansion() = default;

    std::vector<json> hypers{};
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
    std::vector<json> rep_hypers{{{"max_radial", 4}, {"max_angular", 2}},
                                 {{"max_radial", 6}, {"max_angular", 4}}};
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

    SingleHypersSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };

    ~SingleHypersSphericalExpansion() = default;

    std::vector<json> hypers{};
    std::vector<json> fc_hypers{
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.5}, {"unit", "AA"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "AA"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{{{"max_radial", 2},
                                  {"max_angular", 2}}};
  };

  struct SphericalExpansionTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    SphericalExpansionTestData() : Parent{} {
      this->get_ref(this->ref_filename);
    }
    ~SphericalExpansionTestData() = default;
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
  template <typename RepManager>
  class RepresentationManagerGradientCalculator {
   public:
    using Structure_t = AtomicStructure<3>;
    using Key_t = typename RepManager::Key_t;
    using PairRef_t = typename RepManager::Manager_t::template ClusterRef<2>;
    static const size_t n_arguments = 3;

    template <typename T>
    friend class RepresentationManagerGradientFixture;

    RepresentationManagerGradientCalculator(
        RepManager & representation,
        std::shared_ptr<typename RepManager::Manager_t> structure_manager,
        Structure_t atomic_structure)
        : representation{representation}, structure_manager{structure_manager},
          atomic_structure{atomic_structure},
          center_it{structure_manager->begin()} {}

    ~RepresentationManagerGradientCalculator() = default;

    Eigen::Array<double, 1, Eigen::Dynamic>
    f(const Eigen::Ref<const Eigen::Vector3d> & center_position) {
      auto center = *center_it;
      Structure_t modified_structure{this->atomic_structure};
      modified_structure.positions.col(center.get_index()) = center_position;
      this->structure_manager->update(modified_structure);
      representation.compute();
      auto & coeffs_center = representation.expansions_coefficients[center];
      auto keys_center =
          representation.expansions_coefficients.get_keys(center);
      Key_t center_key{center.get_atom_type()};
      size_t n_coeffs_per_key{static_cast<size_t>(
          representation.expansions_coefficients.get_nb_comp())};
      size_t n_coeffs_center{n_coeffs_per_key * keys_center.size()};
      // Packed array containing: The center coefficients (all species) and
      // the neighbour coefficients (only same species as center)
      Eigen::ArrayXd coeffs_pairs(n_coeffs_center +
                                  center.size() * n_coeffs_per_key);

      size_t result_idx{0};
      for (auto & key : keys_center) {
        Eigen::Map<Eigen::RowVectorXd> coeffs_flat(coeffs_center[key].data(),
                                                   n_coeffs_per_key);
        coeffs_pairs.segment(result_idx, n_coeffs_per_key) = coeffs_flat;
        result_idx += n_coeffs_per_key;
      }
      for (auto neigh : center) {
        auto & coeffs_neigh = representation.expansions_coefficients[neigh];
        // The neighbour gradient (i =/= j) only contributes to the channel
        // associated with the _center_ type (the type of the atom that's
        // moving)
        Eigen::Map<Eigen::ArrayXd> coeffs_flat(coeffs_neigh[center_key].data(),
                                               n_coeffs_per_key);
        coeffs_pairs.segment(result_idx, n_coeffs_per_key) = coeffs_flat;
        result_idx += n_coeffs_per_key;
      }

      // Reset the atomic structure for the next iteration
      this->structure_manager->update(this->atomic_structure);
      return coeffs_pairs.transpose();
    }

    Eigen::Array<double, 3, Eigen::Dynamic>
    grad_f(const Eigen::Ref<const Eigen::Vector3d> & /*center_position*/) {
      using Matrix3Xd_RowMaj_t =
          Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
      // Assume f() was already called and updated the position
      // center_it->position() = center_position;
      // representation.compute();
      auto center = *center_it;
      auto keys_center =
          representation.expansions_coefficients.get_keys(center);
      Key_t center_key{center.get_atom_type()};
      size_t n_coeffs_per_key{static_cast<size_t>(
          representation.expansions_coefficients.get_nb_comp())};
      size_t n_coeffs_center{n_coeffs_per_key * keys_center.size()};
      Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>
        grad_coeffs_pairs(3,
                          n_coeffs_center + (center.size() * n_coeffs_per_key));
      auto & grad_coeffs_center =
          representation.expansions_coefficients_gradient[center];
      size_t col_offset{0};
      for (auto & key : keys_center) {
        // Here the 'flattening' retains the 3 Cartesian dimensions as rows,
        // since they vary the slowest within each key
        Eigen::Map<Matrix3Xd_RowMaj_t> grad_coeffs_flat(
            grad_coeffs_center[key].data(), 3, n_coeffs_per_key);
        grad_coeffs_pairs.block(0, col_offset, 3, n_coeffs_per_key) =
            grad_coeffs_flat;
        col_offset += n_coeffs_per_key;
      }
      for (auto neigh : center) {
        typename RepManager::Key_t neigh_key{neigh.get_atom_type()};
        // We need grad_i c^{ji} -- using just 'neigh' would give us
        // grad_j c^{ij}, hence the swap
        // TODO(max,felix) how does this apply to other representations,
        // especially SOAP (sorry, SphericalInvariants<λ=0>)?
        auto neigh_swap{swap_pair_key(neigh)};
        auto & grad_coeffs_neigh =
            representation.expansions_coefficients_gradient[neigh_swap];
        Eigen::Map<Matrix3Xd_RowMaj_t> grad_coeffs_flat(
            grad_coeffs_neigh[center_key].data(), 3, n_coeffs_per_key);
        grad_coeffs_pairs.block(0, col_offset, 3, n_coeffs_per_key) =
            grad_coeffs_flat;
        // The offset keeps advancing from neighbour to neighbour, because the
        // neighbour index has also been flattened out
        col_offset += n_coeffs_per_key;
      }
      return grad_coeffs_pairs;
    }

   private:
    RepManager & representation;
    std::shared_ptr<typename RepManager::Manager_t> structure_manager;
    Structure_t atomic_structure;
    typename RepManager::Manager_t::iterator center_it;

    inline void advance_center() { ++this->center_it; }

    /**
     * Swap a ClusterRef (i, j) so it refers to (j, i) instead
     *
     * @todo wouldn't this be better as a member of StructureManager
     *       (viz. AdaptorNeighbourList<whatever>)?
     */
    PairRef_t swap_pair_key(const PairRef_t & pair_key) {
      // Get the atom index to the corresponding atom tag
      size_t access_index = structure_manager->get_atom_index(pair_key.back());
      auto new_center_it{structure_manager->get_iterator_at(access_index)};
      // Return cluster ref at which the iterator is currently pointing
      auto && new_center{*new_center_it};
      // Iterate until (j,i) is found
      for (auto new_pair : new_center) {
        if (new_pair.back() == pair_key.front()) {
          return new_pair;
        }
      }
      std::stringstream err_str{};
      err_str << "Didn't find symmetric pair for pair (i=" << pair_key.front()
              << ", j=" << pair_key.back() << ").";
      throw std::range_error(err_str.str());
    }
  };

  /**
   * Test fixture holding the gradient calculator and structure manager
   *
   * Holds data (i.e. function values, gradient directions) and iterates through
   * the list of centers
   */
  template <typename RepManager_t>
  class RepresentationManagerGradientFixture : public GradientTestFixture {
   public:
    using StdVector2Dim_t = std::vector<std::vector<double>>;
    using Calculator_t = RepresentationManagerGradientCalculator<RepManager_t>;
    using StructureManager_t = typename RepManager_t::Manager_t;

    static const size_t n_arguments = 3;

    RepresentationManagerGradientFixture(
        std::string filename, std::shared_ptr<StructureManager_t> structure,
        Calculator_t & calc)
        : structure{structure}, center_it{structure->begin()},
          calculator{calc} {
      json input_data;
      std::ifstream input_file{filename};
      input_file >> input_data;

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

    // StdVector2Dim_t function_inputs{};
    // Eigen::MatrixXd displacement_directions{};
    // VerbosityValue verbosity{VerbosityValue::NORMAL};

    std::shared_ptr<StructureManager_t> structure;
    typename StructureManager_t::iterator center_it;
    Calculator_t & calculator;
  };

  struct MultipleStructureSortedCoulomb
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSortedCoulomb() : Parent{} {};
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<json> hypers{{{"central_decay", 0.5},
                              {"interaction_cutoff", 10.},
                              {"interaction_decay", 0.5},
                              {"size", 120},
                              {"sorting_algorithm", "distance"}},
                             {{"central_decay", 0.5},
                              {"interaction_cutoff", 10.},
                              {"interaction_decay", 0.5},
                              {"size", 120},
                              {"sorting_algorithm", "row_norm"}}};
  };

  struct SortedCoulombTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    SortedCoulombTestData() : Parent{} { this->get_ref(this->ref_filename); }
    ~SortedCoulombTestData() = default;

    // name of the file containing the reference data. it has been generated
    // with the following python code:
    // script/generate_sorted_coulomb_ref_data.py

    const bool consider_ghost_neighbours{false};
    std::string ref_filename{"reference_data/sorted_coulomb_reference.ubjson"};
  };

  template <class BaseFixture, template <class> class RepresentationManager>
  struct RepresentationFixture : MultipleStructureFixture<BaseFixture> {
    using Parent = MultipleStructureFixture<BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = RepresentationManager<Manager_t>;

    RepresentationFixture() : Parent{} {}
    ~RepresentationFixture() = default;

    std::list<Representation_t> representations{};
  };

  /* ---------------------------------------------------------------------- */

}  // namespace rascal

#endif  // TESTS_TEST_REPRESENTATION_MANAGER_HH_
