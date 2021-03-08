/**
 * @file   rascal/representations/calculator_spherical_expansion_kspace.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Andrea Grifasi <andrea.grifasi@epfl.ch>
 * @author Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   19 October 2018
 *
 * @brief  Compute the spherical harmonics expansion of the local atom density
 *
 * Copyright © 2018 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_KSPACE_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_KSPACE_HH_

#include "rascal/math/spherical_harmonics.hh"
#include "rascal/math/kvec_generator.hh"
#include "rascal/math/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/cutoff_functions.hh"
#include "rascal/representations/radial_scattering_factors.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/property_block_sparse.hh"
#include "rascal/structure_managers/structure_manager.hh"
#include "rascal/utils/utils.hh"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <memory>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace rascal {

  template <internal::RadialBasisType Type, class Hypers>
  auto make_radial_integral_kspace(const Hypers & basis_hypers) {
    return std::static_pointer_cast<internal::RadialContributionKspaceBase>(
        std::make_shared<internal::RadialContributionKspace<Type>>(basis_hypers));
  }

  template <internal::RadialBasisType Type>
  auto downcast_radial_integral_kspace(
      const std::shared_ptr<internal::RadialContributionKspaceBase> &
          radial_integral) {
    return std::static_pointer_cast<internal::RadialContributionKspace<Type>>(
        radial_integral);
  }

  /**
   * Handles the expansion of an environment in a spherical and radial basis.
   *
   * The local environment of each atom is represented by Gaussians of a
   * certain width (user-defined; can be constant, species-dependent, or
   * radially dependent).  This density field is expanded in an angular basis
   * of spherical harmonics (à la SOAP) and a radial basis of
   * either Gaussians (again, as in SOAP/SphericalInvariants) or one of the more
   * recent bases currently under development.
   */
  class CalculatorKspaceSphericalExpansion : public CalculatorBase {
   public:
    using Parent = CalculatorBase;
    using Hypers_t = typename Parent::Hypers_t;
    using ReferenceHypers_t = Parent::ReferenceHypers_t;
    using Key_t = typename Parent::Key_t;

    template <class StructureManager>
    using Property_t = BlockSparseProperty<double, 1, StructureManager, Key_t>;
    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, StructureManager, Key_t>;

    template <class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;
    template <class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    using Matrix_t = math::Matrix_t;
    using Vector_t = math::Vector_t;
    using Matrix_Ref = math::Matrix_Ref;
    using Vector_Ref = math::Vector_Ref;

    /**
     * Set the hyperparameters of this descriptor from a json-like container.
     *
     * @param hypers structure (usually parsed from json) containing the
     *               options and hyperparameters
     *
     * @todo (max, felix) document the SOAP/SphExpn-specific hypers here
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the structure
     */
    void set_hyperparameters(const Hypers_t & hypers) override {
      using internal::RadialBasisType;
      this->hypers = hypers;

      this->max_radial = hypers.at("max_radial").get<size_t>();
      this->max_angular = hypers.at("max_angular").get<size_t>();

      auto smearing_hypers = hypers.at("gaussian_density").get<json>();
      this->smearing = smearing_hypers.at("gaussian_sigma").at("value").get<double>();

      if (hypers.count("compute_gradients")) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }

      if (hypers.count("expansion_by_species_method")) {
        std::set<std::string> possible_expansion_by_species{
            {"environment wise", "user defined", "structure wise"}};
        auto expansion_by_species_tmp =
            hypers.at("expansion_by_species_method").get<std::string>();
        if (possible_expansion_by_species.count(expansion_by_species_tmp)) {
          this->expansion_by_species = expansion_by_species_tmp;
        } else {
          std::stringstream err_str{};
          err_str << "expansion_by_species_method provided:'"
                  << expansion_by_species_tmp
                  << "' is not part of the implemented methods: '";
          for (const auto & val : possible_expansion_by_species) {
            err_str << val << "', ";
          }
          throw std::logic_error(err_str.str());
        }
      } else {
        // default value for backward compatibility
        this->expansion_by_species = "environment wise";
      }

      if (hypers.count("global_species")) {
        auto species = hypers.at("global_species").get<Key_t>();
        for (const auto & sp : species) {
          this->global_species.insert({sp});
        }
      } else {
        if (this->expansion_by_species == "user defined") {
          std::stringstream err_str{};
          err_str << "expansion_by_species_method is 'user defined'"
                  << " but global_species is not defined.";
          throw std::logic_error(err_str.str());
        }
        this->global_species.clear();
      }

      this->spherical_harmonics.precompute(this->max_angular,this->compute_gradients);

      // create the class that will compute the radial terms of the
      // expansion. the atomic smearing is an integral part of the
      // radial contribution
      auto radial_contribution_hypers =
          hypers.at("radial_contribution").get<json>();
      auto radial_contribution_type =
          radial_contribution_hypers.at("type").get<std::string>();

      // create the class that will compute the radial terms of the
      // expansion. the atomic smearing is an integral part of the
      // radial contribution
      if (radial_contribution_type == "GTO") {
        auto rc_shared = std::make_shared<
            internal::RadialContributionKspace<RadialBasisType::GTO>>(hypers);
        this->radial_integral = rc_shared;
        this->radial_integral_type = RadialBasisType::GTO;
      } else {
        throw std::logic_error("Requested Radial contribution type \'" +
                               radial_contribution_type +
                               "\' has not been implemented.  Must be one of" +
                               ": \'GTO\'. ");
      }

      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      // auto fc_type = fc_hypers.at("type").get<std::string>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
      //this->cutoff_smooth_width = fc_hypers.at("smooth_width").at("value");
      this->set_name(hypers);
    }

    bool operator==(const CalculatorKspaceSphericalExpansion & other) const {
      bool is_equal{
          (this->does_gradients() == other.does_gradients()) and
          (this->interaction_cutoff == other.interaction_cutoff) and
          (this->interpolator_accuracy == other.interpolator_accuracy) and
          (this->max_radial == other.max_radial) and
          (this->max_angular == other.max_angular) and
          (this->radial_integral_type == other.radial_integral_type) and
          (this->smearing == other.smearing)};
      return is_equal;
    }

    /**
     * Returns if the calculator is able to compute gradients of the
     * representation w.r.t. atomic positions ?
     */
    bool does_gradients() const override { return this->compute_gradients; }

    /**
     * Construct a new Calculator using a hyperparameters container
     *
     * See set_hyperparameters() for a description of the hypers
     *
     * @todo (max) ffs, why isn't the link above working in Sphinx?
     *             And why aren't the todos showing?!
     *
     * @param hyper container (usually parsed from json) for the options and
     *              hyperparameters
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the container
     */
    explicit CalculatorKspaceSphericalExpansion(const Hypers_t & hyper)
        : CalculatorBase{} {
      this->set_default_prefix("spherical_expansion_");
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    CalculatorKspaceSphericalExpansion(const CalculatorKspaceSphericalExpansion & other) =
        delete;

    //! Move constructor
    CalculatorKspaceSphericalExpansion(CalculatorKspaceSphericalExpansion && other) noexcept
        : CalculatorBase{std::move(other)}, interaction_cutoff{std::move(
                                                other.interaction_cutoff)},
          smearing{std::move(other.smearing)},
          max_radial{std::move(other.max_radial)}, max_angular{std::move(
                                                       other.max_angular)},
          compute_gradients{std::move(other.compute_gradients)},
          expansion_by_species{std::move(other.expansion_by_species)},
          global_species{std::move(other.global_species)},
          radial_integral{std::move(other.radial_integral)},
          radial_integral_type{std::move(other.radial_integral_type)},
          spherical_harmonics{std::move(other.spherical_harmonics)} {}

    //! Destructor
    virtual ~CalculatorKspaceSphericalExpansion() = default;

    //! Copy assignment operator
    CalculatorKspaceSphericalExpansion &
    operator=(const CalculatorKspaceSphericalExpansion & other) = delete;

    //! Move assignment operator
    CalculatorKspaceSphericalExpansion &
    operator=(CalculatorKspaceSphericalExpansion && other) = default;

    /**
     * Compute representation for a given structure manager.
     *
     * @tparam StructureManager a (single or collection)
     * of structure manager(s) (in an iterator) held in shared_ptr
     */
    template <class StructureManager>
    void compute(StructureManager & managers);

    /**
     * loop over a collection of manangers if it is an iterator.
     * Or just call compute_impl() if it's a single manager (see below)
     */
    template <
        internal::RadialBasisType RadialType,
        class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<RadialType>(manager);
      }
    }

    //! single manager case
    template <internal::RadialBasisType RadialType,
              class StructureManager,
              std::enable_if_t<
                  not(internal::is_proper_iterator<StructureManager>::value),
                  int> = 0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<RadialType>(manager);
    }

    //! Compute the spherical exansion given several options
    template <internal::RadialBasisType RadialType,
              class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

   protected:
    //! cutoff radius r_c defining the radius around the atom for expanding the potential
    double interaction_cutoff{};
    //! defines the maximal mean error allowed to the interpolator when fitting
    //! the reference
    double interpolator_accuracy{};
    //! smearing
    double smearing{};
    //! number of radial basis function to use in the expansion
    size_t max_radial{};
    /**
     * number of angular channels used in the expansion, i.e. all Y_l^m with
     * l < max_angular + 1 are used (the +1 is to follow GAP's convention).
     */
    size_t max_angular{};
    //! controls the computation of the gradients of the expansion wrt. atomic
    //! positions
    bool compute_gradients{};
    /**
     * defines the method to determine the set of species to use in the
     * expansion
     */
    std::string expansion_by_species{};

    //! user defined species appearing in the expansion indexing
    std::set<Key_t> global_species{};

    std::shared_ptr<internal::RadialContributionKspaceBase> radial_integral{};
    internal::RadialBasisType radial_integral_type{};

    math::SphericalHarmonics spherical_harmonics{};

    /**
     * set up chemical keys of the expension so that only species appearing in
     * the environment are present and initialize coeffs to zero.
     *
     * For gradients associated with pair_ii the keys will be the ones in the
     * environment but the ones associated with pair_ij will only contain the
     * non zero keys.
     */
    template <class StructureManager>
    void initialize_expansion_environment_wise(
        std::shared_ptr<StructureManager> & managers,
        Property_t<StructureManager> & expansions_coefficients,
        PropertyGradient_t<StructureManager> &
            expansions_coefficients_gradient);

    /**
     * set up chemical keys of the expension so that all species in the
     * structure will be used in the as keys for the expansion and initialize
     * coeffs to zero.
     * For gradients associated with pair_ii the keys will be the one of the
     * structure but the ones associated with pair_ij will only contain the
     * non zero keys.
     *
     * @throw runtime_error when all the species of the structure are not
     * present in global_species
     */
    template <class StructureManager>
    void initialize_expansion_structure_wise(
        std::shared_ptr<StructureManager> & managers,
        Property_t<StructureManager> & expansions_coefficients,
        PropertyGradient_t<StructureManager> &
            expansions_coefficients_gradient);

    /**
     * Set up chemical keys of the expension using global_species for the keys
     * appearing in the expansion and initialize coeffs to zero.
     * For gradients associated with pair_ii the keys will be the one of
     * global_species but the ones associated with pair_ij will only contain the
     * non zero keys.
     *
     * @throw runtime_error when all the species of the structure are not
     * present in global_species
     */
    template <class StructureManager>
    void initialize_expansion_with_global_species(
        std::shared_ptr<StructureManager> & managers,
        Property_t<StructureManager> & expansions_coefficients,
        PropertyGradient_t<StructureManager> &
            expansions_coefficients_gradient);
  };

  // compute classes template construction
  template <class StructureManager>
  void CalculatorKspaceSphericalExpansion::compute(StructureManager & managers) {
    // We don't specialize on cutoff function anymore, just radial contribution
    // (and not even that, really, since we only implement GTO)
    using internal::RadialBasisType;
    switch(this->radial_integral_type) {
        case RadialBasisType::GTO: {
            this->compute_loop<RadialBasisType::GTO>(managers);
            break;
        }
        default:
          // The control flow really should never reach here...
          std::basic_ostringstream<char> err_message;
          err_message << "Invalid radial basis type encountered ";
          err_message << "(This is a bug.  Debug info for developers: "
                      << "radial_integral_type == ";
          err_message << static_cast<int>(this->radial_integral_type);
          err_message << ")" << std::endl;
          throw std::logic_error(err_message.str());
        }
  }

  /**
   * Compute the spherical expansion
   */
  template <internal::RadialBasisType RadialType,
            class StructureManager>
  void CalculatorKspaceSphericalExpansion::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using Prop_t = Property_t<StructureManager>;
    using PropGrad_t = PropertyGradient_t<StructureManager>;
    constexpr static bool IsHalfNL{
        StructureManager::traits::NeighbourListType ==
        AdaptorTraits::NeighbourListType::half};
    using math::PI;
    using math::pow;
    constexpr bool ExcludeGhosts{true};
    const bool is_not_masked{manager->is_not_masked()};
    const bool compute_gradients{this->compute_gradients};
    if (not is_not_masked and compute_gradients) {
      throw std::logic_error("Can't compute spherical expansion gradients with "
                             "masked center atoms");
    }
    if (not is_not_masked and IsHalfNL) {
      std::stringstream err_str{};
      err_str << "Half neighbor list should only be used when all the "
              << "atoms inside the unit cell are centers, i.e. "
              << "center_atoms_mask should not mask atoms.";
      throw std::runtime_error(err_str.str());
    }
    auto manager_root = extract_underlying_manager<0>(manager);
    auto cell_length = manager_root->get_cell_length();
    auto cell = manager_root->get_cell();
    auto tcoords = manager_root->get_positions();
    Matrix_t coords = tcoords.transpose();
    int natoms = coords.rows();
    const double volume = cell.determinant();

    // Radius of the sphere in reciprocal space defining the maximum spatial resolution
    // pi/sigma was shown to be enough to converge the density field in TENSOAP
    double kcut = 2.0*PI/(this->smearing);
    
    // get cell vectors in reciprocal space
    Matrix_t tcell = cell.transpose();
    Matrix_t bvecs = 2.0*PI * tcell.inverse();

    /* PREPARATION:
    Determine the optimal bounds n1max, n2max, n3max that define the search space box
    Roughly speaking, our goal will then be to find all vectors k = n1*b1 + n2*b2 + n3*b3,
    where 0<n1<=n1max, abs(n2)<=n2max, abs(n3)<=n3max with some special care at the boundaries.
    The derivation of the formulae is found in a supporting document.
    */

    // Compute the matrix M_ij = b_i * b_j, the representation matrix of the inner product in index space.
    double M11 = bvecs(0,0) * bvecs(0,0) + bvecs(0,1) * bvecs(0,1) + bvecs(0,2) * bvecs(0,2);
    double M22 = bvecs(1,0) * bvecs(1,0) + bvecs(1,1) * bvecs(1,1) + bvecs(1,2) * bvecs(1,2);
    double M33 = bvecs(2,0) * bvecs(2,0) + bvecs(2,1) * bvecs(2,1) + bvecs(2,2) * bvecs(2,2);
    double M12 = bvecs(0,0) * bvecs(1,0) + bvecs(0,1) * bvecs(1,1) + bvecs(0,2) * bvecs(1,2);
    double M13 = bvecs(0,0) * bvecs(2,0) + bvecs(0,1) * bvecs(2,1) + bvecs(0,2) * bvecs(2,2);
    double M23 = bvecs(1,0) * bvecs(2,0) + bvecs(1,1) * bvecs(2,1) + bvecs(1,2) * bvecs(2,2);

    // Get the optimal boundaries for the search space box (formula from supporting information)
    double Mbar = M11 * M22 * M33 - (M11 * M23 * M23 + M22 * M13 * M13 + M33 * M12 * M12) + 2 * M12 * M13 * M23;
    int n1max = floor(sqrt((M22 * M33 - M23 * M23) / Mbar) * kcut);
    int n2max = floor(sqrt((M11 * M33 - M13 * M13) / Mbar) * kcut);
    int n3max = floor(sqrt((M11 * M22 - M12 * M12) / Mbar) * kcut);

    // Total number of points that will be checked
    int numtot = n3max + n2max * (2 * n3max + 1) + n1max * (2 * n2max + 1) * (2 * n3max + 1);

    // Initialize k-vectors generation
    math::Kvectors k_vectors(numtot);

    // Compute k-vectors on a semisphere in the space spanned by n1,n2,n3 
    k_vectors.precompute(n1max,n2max,n3max,bvecs,kcut);
    int n_k = k_vectors.nk; // number of k-vectors
    Vector_t k_val = k_vectors.kval; // k-moduli
    Matrix_t k_vec = k_vectors.kvec; // k-vectors

    // Console output to check code (TODO: DELETE THIS AFTER TESTING PHASE)
    double succratio = ((double)(n_k)) / numtot;
    std::cout << "Number of k-vectors = " << n_k << "\n";
    std::cout << "Total number of points = " << numtot << "\n";
    std::cout << "Ratio of successful points = " << succratio << "\n";
    std::cout << "Ideal success ratio of circle = " << 3.1415 / 6 << "\n";

    auto pbc = manager_root->get_periodic_boundary_conditions();
    bool is_cutoff_too_large{false};
    for (size_t i_dim{0}; i_dim < ThreeD; ++i_dim) {
      if (pbc[i_dim]) {
        if (cell_length[i_dim] < 2. * this->interaction_cutoff) {
          is_cutoff_too_large = true;
        }
      }
    }
    if (IsHalfNL and is_cutoff_too_large) {
      std::stringstream err_str{};
      err_str << "Half neighbor list should only be used when the diameter of "
              << "the spherical expansion is smaller than the unit cell "
              << "in periodic directions: "
              << "[" << cell_length.transpose() << "] > "
              << 2 * this->interaction_cutoff;
      throw std::runtime_error(err_str.str());
    }

    auto && expansions_coefficients{*manager->template get_property<Prop_t>(
        this->get_name(), true, true, ExcludeGhosts)};

    auto && expansions_coefficients_gradient{
        *manager->template get_property<PropGrad_t>(this->get_gradient_name(),
                                                    true, true)};

    // if the representation has already been computed for the current
    // structure then do nothing
    if (expansions_coefficients.is_updated()) {
      return;
    }

    // downcast cutoff and radial contributions so they are functional
    auto radial_integral{
        downcast_radial_integral_kspace<RadialType>(
            this->radial_integral)};

    auto n_row{this->max_radial};
    // to store linearly all l,m components with
    // -l-1<=m<=l+1 needs (l+1)**2 elements
    auto n_col{(this->max_angular + 1) * (this->max_angular + 1)};
    expansions_coefficients.clear();
    expansions_coefficients.set_shape(n_row, n_col);

    if (compute_gradients) {
      expansions_coefficients_gradient.clear();
      // Row-major ordering, so the Cartesian (spatial) index varies slowest
      expansions_coefficients_gradient.set_shape(ThreeD * n_row, n_col);
    }

    if (this->expansion_by_species == "environment wise") {
      this->initialize_expansion_environment_wise(
          manager, expansions_coefficients, expansions_coefficients_gradient);
    } else if (this->expansion_by_species == "user defined") {
      this->initialize_expansion_with_global_species(
          manager, expansions_coefficients, expansions_coefficients_gradient);
    } else if (this->expansion_by_species == "structure wise") {
      this->initialize_expansion_structure_wise(
          manager, expansions_coefficients, expansions_coefficients_gradient);
    } else {
      throw std::runtime_error("should not arrive here");
    }

    // Fourier components of a Gaussian charge 
    auto G_k = math::Vector_t(n_k); 
    
    // Spherical harmonics matrix N_k * (l_max+1)^2
    auto Y_lm = math::Matrix_t(n_k, n_col); 
    
    // Radial integrals N_k * (n_max*(l_max+1))
    auto I_nl = math::Matrix_t(n_k, this->max_radial*(this->max_angular+1)); 
    
    // Trigonometric matrices N_k * N_atoms 
    auto cos_ki = math::Matrix_t(n_k,natoms); 
    auto sin_ki = math::Matrix_t(n_k,natoms); 
    
    // precompute quantities needed 
    for (int ik{0}; ik < n_k; ++ik) {
   
      // Fourier charge at k_val
      G_k(ik) = std::exp(-0.5 * pow(this->smearing*k_val(ik) , 2) ) ;/// pow(k_val(ik),2);

      // Harmonics at k_dir
      if (k_val(ik)!=0.0){
        Eigen::Vector3d k_dir;
        k_dir(0) = k_vec(ik,0)/k_val(ik);
        k_dir(1) = k_vec(ik,1)/k_val(ik);
        k_dir(2) = k_vec(ik,2)/k_val(ik);
        this->spherical_harmonics.calc(k_dir);
        auto && harmonics{spherical_harmonics.get_harmonics()};
        for (size_t lm{0}; lm < n_col; ++lm) {
          Y_lm(ik,lm) = harmonics(lm);
        }
      } else {
	Y_lm(ik,0) = 0.5/std::sqrt(4*PI); 
      }

      // Radial integral at k_val
      auto && radint = radial_integral->compute_radial_integral(k_val(ik)); // n_max * (l_max+1)
      int nl_size = this->max_radial*(this->max_angular+1); 
      for (int nl{0}; nl < nl_size; ++nl){
        I_nl(ik,nl) = radint(nl);
      }

      // Trigonometric functions
      for (int iat{0}; iat < natoms; ++iat) {
	double arg = coords(iat,0) * k_vec(ik,0) 
		   + coords(iat,1) * k_vec(ik,1) 
		   + coords(iat,2) * k_vec(ik,2);
        cos_ki(ik,iat) = std::cos(arg); 
        sin_ki(ik,iat) = std::sin(arg); 
      } 

    }
    
    // Start the accumulation
    for (auto center : manager) {

      auto atom_i_tag = center.get_atom_tag();
      size_t iat = manager->get_atom_index(atom_i_tag); 
      // c^{i}
      auto & coefficients_center = expansions_coefficients[center];
      // \grad_i c^{i}
      //auto & coefficients_center_gradient =
      //    expansions_coefficients_gradient[center.get_atom_ii()];
      Key_t center_type{center.get_atom_type()};

      // loop over k-vectors
      for (int ik{0}; ik < n_k; ++ik) {
        
        // Start the accumulation with the central atom contribution
        int nl_idx{0};
        for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          int l_block_idx{0};
          for (size_t angular_l{0}; angular_l < this->max_angular+1; ++angular_l) {
            // l odd contributions vanish by k-symmetry for central atom
	    int size_m = 2*angular_l+1;
            if (angular_l%2==0) {
	      for (int mval{0}; mval < size_m; ++mval) {
	        coefficients_center[center_type](radial_n,l_block_idx+mval) += 
	          I_nl(ik,nl_idx) * Y_lm(ik,l_block_idx+mval) * G_k(ik) 
		  * 16.0 * pow(PI,2) / volume; 
	      }
            }
            l_block_idx += size_m;
            nl_idx += 1;
          }
	}

        // loop over atoms over the entire system using centers idx once again      
        for (auto neigh : manager) {
          
	  auto atom_j_tag = neigh.get_atom_tag();
          size_t jat = manager->get_atom_index(atom_j_tag);
	  // accumulate for atoms different from the central
	  if (jat != iat) {
            Key_t neigh_type{neigh.get_atom_type()};

            // compute the coefficients
            int nl_idx = 0;
            for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
              int l_block_idx = 0;
              for (size_t angular_l{0}; angular_l < this->max_angular+1; ++angular_l) {
	        double phase_factor;
                if (angular_l%2==0) {
	         // real phase factor for l even
	         phase_factor = pow(-1.0,angular_l/2) 
	             * (cos_ki(ik,jat)*cos_ki(ik,iat) + sin_ki(ik,jat)*sin_ki(ik,iat));
                } else {
	         // real phase factor for l odd 
	         phase_factor = -pow(-1.0,(angular_l+1)/2) 
	             * (sin_ki(ik,jat)*cos_ki(ik,iat) - cos_ki(ik,jat)*sin_ki(ik,iat));
	        } 
	        int size_m = 2*angular_l+1;
	        for (int mval{0}; mval < size_m; ++mval) {
	          coefficients_center[neigh_type](radial_n,l_block_idx+mval) += 
                    I_nl(ik,nl_idx) * Y_lm(ik,l_block_idx+mval) * G_k(ik)
                    * phase_factor * 16.0 * pow(PI,2) / volume;	
                }
                l_block_idx += size_m;
                nl_idx += 1;
              }
	    }
	  }

//          // half list branch for c^{ji} terms using
//          // c^{ij}_{nlm} = (-1)^l c^{ji}_{nlm}.
//          if (IsHalfNL) {
//            if (is_center_atom) {
//              auto & coefficients_neigh{expansions_coefficients[atom_j]};
//              auto coefficients_neigh_by_type{coefficients_neigh[center_type]};
//              int l_block_idx = 0;
//              double parity{1.};
//              for (size_t angular_l{0}; angular_l < this->max_angular + 1;
//                   ++angular_l) {
//                size_t l_block_size{2 * angular_l + 1};
//                coefficients_neigh_by_type.block(0, l_block_idx, max_radial,
//                                                 l_block_size) +=
//                    parity *
//                    c_ij_nlm.block(0, l_block_idx, max_radial, l_block_size);
//                l_block_idx += l_block_size;
//                parity *= -1.;
//              }
//            }
//          }

//          // compute the gradients of the coefficients with respect to atoms positions
//          // but only if the neighbour is not an image of the center!
//          // (the periodic images move with the center, so their contribution to the center gradient is zero)
//          if (compute_gradients and (atom_j_tag != atom_i_tag)) {  // NOLINT
//            // \grad_j c^i
//            auto & coefficients_neigh_gradient =
//                expansions_coefficients_gradient[neigh];
//
//            // The type of the contribution c^{ij} to the coefficient c^{i}
//            // depends on the type of j (and it is the same for the gradients)
//            // In the following atom i is of type a and atom j is of type b
//
//            // grad_i c^{ib}
//            auto && gradient_center_by_type{
//                coefficients_center_gradient[neigh_type]};
//            // grad_j c^{ib}
//            auto && gradient_neigh_by_type{
//                coefficients_neigh_gradient[neigh_type]};
//
//            // grad_j c^{ij}
//            Matrix_t pair_gradient_contribution{this->max_radial,
//                                                this->max_angular + 1};
//            for (int cartesian_idx{0}; cartesian_idx < ThreeD;
//                   ++cartesian_idx) {
//              l_block_idx = 0;
//              for (size_t angular_l{0}; angular_l < this->max_angular + 1;
//                  ++angular_l) {
//                size_t l_block_size{2 * angular_l + 1};
//                pair_gradient_contribution.resize(this->max_radial, l_block_size);
//                // grad_i c^{ib} = - \sum_{j} grad_j c^{ijb}
//                gradient_center_by_type.block(
//                    cartesian_idx * max_radial, l_block_idx,
//                    max_radial, l_block_size) -= pair_gradient_contribution;
//                // grad_j c^{ib} =  grad_j c^{ijb}
//                gradient_neigh_by_type.block(
//                    cartesian_idx * max_radial, l_block_idx,
//                    max_radial, l_block_size) = pair_gradient_contribution;
//                l_block_idx += l_block_size;
//                // clang-format on
//              }  // for (angular_l)
//            }    // for cartesian_idx
//
//            // half list branch for accumulating parts of grad_j c^{j} using
//            // grad_j c^{ji a} = (-1)^l grad_j c^{ij b}
//            if (IsHalfNL) {
//              if (is_center_atom) {
//                // grad_j c^{j}
//                auto & coefficients_neigh_center_gradient =
//                    expansions_coefficients_gradient[neigh.get_atom_jj()];
//                // grad_j c^{j a}
//                auto gradient_neigh_center_by_type =
//                    coefficients_neigh_center_gradient[center_type];
//
//                for (int cartesian_idx{0}; cartesian_idx < ThreeD;
//                     ++cartesian_idx) {
//                  l_block_idx = 0;
//                  double parity{1};
//                  for (size_t angular_l{0}; angular_l < this->max_angular + 1;
//                       ++angular_l) {
//                    size_t l_block_size{2 * angular_l + 1};
//                    // clang-format off
//                    gradient_neigh_center_by_type.block(
//                        cartesian_idx * max_radial, l_block_idx,
//                        max_radial, l_block_size) += parity *
//                                    gradient_neigh_by_type.block(
//                                      cartesian_idx * max_radial, l_block_idx,
//                                      max_radial, l_block_size);
//
//                    l_block_idx += l_block_size;
//                    parity *= -1.;
//                    // clang-format on
//                  }  // for (angular_l)
//                }    // for cartesian_idx
//              }      // if (is_center_atom)
//            }        // if (IsHalfNL)
//          }          // if (compute_gradients)
        }            // for (neigh : center)
      
      }            // for (k-vectors)

      // Normalize and orthogonalize the radial coefficients
      radial_integral->finalize_coefficients(coefficients_center);
      //if (compute_gradients) {
      //  radial_integral->template finalize_coefficients_der<ThreeD>(
      //      expansions_coefficients_gradient, center);
      // }
    }  // for (center : manager)
  }    // compute()


  // STRUCTURE MANAGER STUFF BELOW
  template <class StructureManager>
  void CalculatorKspaceSphericalExpansion::initialize_expansion_environment_wise(
      std::shared_ptr<StructureManager> & manager,
      Property_t<StructureManager> & expansions_coefficients,
      PropertyGradient_t<StructureManager> & expansions_coefficients_gradient) {
    constexpr static bool IsHalfNL{
        StructureManager::traits::NeighbourListType ==
        AdaptorTraits::NeighbourListType::half};
    std::vector<std::set<Key_t>> keys_list{};
    std::vector<std::set<Key_t>> keys_list_grad{};
    std::map<int, int> center_tag2idx{};
    const bool compute_gradients{this->compute_gradients};
    int i_center{0};
    for (auto center : manager) {
      center_tag2idx[center.get_atom_tag()] = i_center;
      i_center++;
      keys_list.emplace_back();
      if (compute_gradients) {
        for (auto neigh : center.pairs_with_self_pair()) {
          (void)neigh;  // to avoid compiler warning
          keys_list_grad.emplace_back();
        }
      }
    }
    int i_grad{0};
    i_center = 0;
    for (auto center : manager) {
      Key_t center_type{center.get_atom_type()};
      auto atom_i_tag = center.get_atom_tag();

      for (auto neigh : center.pairs()) {
        keys_list[i_center].insert({neigh.get_atom_type()});
        if (manager->is_center_atom(neigh) and IsHalfNL) {
          auto atom_j = neigh.get_atom_j();
          auto j_center = center_tag2idx[atom_j.get_atom_tag()];
          keys_list[j_center].insert(center_type);
        }
      }
      keys_list[i_center].insert({center_type});
      if (compute_gradients) {
        keys_list_grad[i_grad].insert(keys_list[i_center].begin(),
                                      keys_list[i_center].end());
        i_grad++;
        for (auto neigh : center.pairs()) {
          auto && atom_j = neigh.get_atom_j();
          auto atom_j_tag = atom_j.get_atom_tag();
          if (atom_j_tag != atom_i_tag) {
            Key_t neigh_type{neigh.get_atom_type()};
            keys_list_grad[i_grad].insert(neigh_type);
          }
          i_grad++;
        }
      }  // if (compute_gradients)
      i_center++;
    }

    expansions_coefficients.resize(keys_list);
    expansions_coefficients.setZero();

    if (compute_gradients) {
      expansions_coefficients_gradient.resize(keys_list_grad);
      expansions_coefficients_gradient.setZero();
    } else {
      expansions_coefficients_gradient.resize();
    }
  }

  template <class StructureManager>
  void CalculatorKspaceSphericalExpansion::initialize_expansion_structure_wise(
      std::shared_ptr<StructureManager> & manager,
      Property_t<StructureManager> & expansions_coefficients,
      PropertyGradient_t<StructureManager> & expansions_coefficients_gradient) {
    std::set<Key_t> keys{};
    for (auto center : manager) {
      Key_t center_type{center.get_atom_type()};
      keys.insert({center_type});
      // there might be masked atoms having different types from the centers
      // so also need to loop over the pairs here
      for (auto neigh : center.pairs()) {
        keys.insert({neigh.get_atom_type()});
      }
    }

    std::vector<std::set<Key_t>> keys_list{};
    std::vector<std::set<Key_t>> keys_list_grad{};
    for (auto center : manager) {
      // Key_t center_type{center.get_atom_type()};
      auto atom_i_tag = center.get_atom_tag();
      keys_list.emplace_back(keys);
      if (this->compute_gradients) {
        keys_list_grad.emplace_back(keys);
        for (auto neigh : center.pairs()) {
          auto && atom_j = neigh.get_atom_j();
          auto atom_j_tag = atom_j.get_atom_tag();
          std::set<Key_t> neigh_types{};
          if (atom_j_tag != atom_i_tag) {
            Key_t neigh_type{neigh.get_atom_type()};
            neigh_types.insert(neigh_type);
          }
          keys_list_grad.emplace_back(neigh_types);
        }
      }
    }

    expansions_coefficients.resize(keys_list);
    expansions_coefficients.setZero();

    if (this->compute_gradients) {
      expansions_coefficients_gradient.resize(keys_list_grad);
      expansions_coefficients_gradient.setZero();
    } else {
      expansions_coefficients_gradient.resize();
    }
  }

  template <class StructureManager>
  void CalculatorKspaceSphericalExpansion::initialize_expansion_with_global_species(
      std::shared_ptr<StructureManager> & manager,
      Property_t<StructureManager> & expansions_coefficients,
      PropertyGradient_t<StructureManager> & expansions_coefficients_gradient) {
    std::vector<std::set<Key_t>> keys_list{};
    std::vector<std::set<Key_t>> keys_list_grad{};

    // check that all species in the structure are present in global_species
    Key_t keys{};
    for (auto center : manager) {
      typename Key_t::value_type center_type{center.get_atom_type()};
      keys.push_back(center_type);
    }
    Key_t missing_keys{};
    for (const auto & key : keys) {
      if (not internal::is_element_in(key, this->global_species)) {
        missing_keys.emplace_back(key);
      }
    }
    if (missing_keys.size() > 0) {
      std::stringstream err_str{};
      err_str << "global_species is missing at least these species: '";
      for (const auto & key : missing_keys) {
        err_str << key << ", ";
      }
      err_str << "'.";
      throw std::runtime_error(err_str.str());
    }

    // build the species list
    for (auto center : manager) {
      // Key_t center_type{center.get_atom_type()};
      keys_list.emplace_back(this->global_species);
      if (this->compute_gradients) {
        keys_list_grad.emplace_back(this->global_species);
        auto atom_i_tag = center.get_atom_tag();
        for (auto neigh : center.pairs()) {
          auto && atom_j = neigh.get_atom_j();
          auto atom_j_tag = atom_j.get_atom_tag();
          std::set<Key_t> neigh_types{};
          if (atom_j_tag != atom_i_tag) {
            Key_t neigh_type{neigh.get_atom_type()};
            neigh_types.insert(neigh_type);
          }
          keys_list_grad.emplace_back(neigh_types);
        }
      }
    }

    expansions_coefficients.resize(keys_list);
    expansions_coefficients.setZero();

    if (this->compute_gradients) {
      expansions_coefficients_gradient.resize(keys_list_grad);
      expansions_coefficients_gradient.setZero();
    } else {
      expansions_coefficients_gradient.resize();
    }
  }

}  // namespace rascal

namespace nlohmann {
  /**
   * Special specialization of the json serialization for non default
   * constructible type.
   */
  template <>
  struct adl_serializer<rascal::CalculatorKspaceSphericalExpansion> {
    static rascal::CalculatorKspaceSphericalExpansion from_json(const json & j) {
      return rascal::CalculatorKspaceSphericalExpansion{j};
    }

    static void to_json(json & j,
                        const rascal::CalculatorKspaceSphericalExpansion & t) {
      j = t.hypers;
    }
  };
}  // namespace nlohmann

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_KSPACE_HH_
