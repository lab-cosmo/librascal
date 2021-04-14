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

#include "rascal/math/kvec_generator.hh"
#include "rascal/math/spherical_harmonics.hh"
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
#include <fstream>

namespace rascal {

  template <internal::RadialBasisType Type, class Hypers>
  auto make_radial_integral_kspace(const Hypers & basis_hypers) {
    return std::static_pointer_cast<internal::RadialContributionKspaceBase>(
        std::make_shared<internal::RadialContributionKspace<Type>>(
            basis_hypers));
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
      this->smearing =
          smearing_hypers.at("gaussian_sigma").at("value").get<double>();

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

      // @memo: test kspace variant that never needs derivatives of
      //        spherical harmonics
      this->spherical_harmonics.precompute(this->max_angular,false);
      //this->spherical_harmonics.precompute(this->max_angular,this->compute_gradients);

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
        throw std::logic_error(
            "Requested Radial contribution type \'" + radial_contribution_type +
            "\' has not been implemented.  Must be one of" + ": \'GTO\'. ");
      }

      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      // auto fc_type = fc_hypers.at("type").get<std::string>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
      // this->cutoff_smooth_width = fc_hypers.at("smooth_width").at("value");
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
    CalculatorKspaceSphericalExpansion(
        const CalculatorKspaceSphericalExpansion & other) = delete;

    //! Move constructor
    CalculatorKspaceSphericalExpansion(
        CalculatorKspaceSphericalExpansion && other) noexcept
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
        internal::RadialBasisType RadialType, class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<RadialType>(manager);
      }
    }

    //! single manager case
    template <internal::RadialBasisType RadialType, class StructureManager,
              std::enable_if_t<
                  not(internal::is_proper_iterator<StructureManager>::value),
                  int> = 0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<RadialType>(manager);
    }

    //! Compute the spherical exansion given several options
    template <internal::RadialBasisType RadialType, class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

   protected:
    //! cutoff radius r_c defining the radius around the atom for expanding the
    //! potential
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
  void
  CalculatorKspaceSphericalExpansion::compute(StructureManager & managers) {
    // We don't specialize on cutoff function anymore, just radial contribution
    // (and not even that, really, since we only implement GTO)
    using internal::RadialBasisType;
    switch (this->radial_integral_type) {
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
  template <internal::RadialBasisType RadialType, class StructureManager>
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

    /*
      -------------------------------------------------------------------
      -------------------------------------------------------------------
      ----------START STEPS DEPENDENT ON CELL BUT NOT ON POSITIONS ------
      -------------------------------------------------------------------
      -------------------------------------------------------------------
    */

    // Define cell and some of the cell dependent parameters
    auto manager_root = extract_underlying_manager<0>(manager);
    auto cell_length = manager_root->get_cell_length();
    auto cell = manager_root->get_cell();
    const double volume = cell.determinant();

    /*
            ----------------------------------------------------------------------------------
            DO SOME PREPARATION ???:

            Does a lot of small checks and auxiliary definitions.
            For some of the steps, not exactly sure what they are doing.
            ----------------------------------------------------------------------------------

    */

    // Initialize some of the required parameters
    // @memo: since we always work with periodic boundary conditions
    //        for the kspace implementation, this could be
    //        replaced by some assertion e.g. assert(pbc=true)
    //        Also, does the restriction regarding the
    //        interaction_cutoff still apply?
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

    // Initialize expansion coefficients and their gradients 
    auto && expansions_coefficients{*manager->template get_property<Prop_t>(
        this->get_name(), true, true, ExcludeGhosts)};

    auto && expansions_coefficients_gradient{
        *manager->template get_property<PropGrad_t>(this->get_gradient_name(),
                                                    true, true)};

    // if the representation has already been computed for the current
    // structure then do nothing
    // @memo: shouldn't this come a little earlier? 
    if (expansions_coefficients.is_updated()) {
      return;
    }

    // downcast cutoff and radial contributions so they are functional
    auto radial_integral{
        downcast_radial_integral_kspace<RadialType>(this->radial_integral)};

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

    /*
      -------------------------------------------------------------------
      GENERATE K-VECTORS GRID:

      Generate the k-space vectors required for the (discrete)
      Fourier transform using the Kvectors class (see ../math)
      The class requires:
      - a cutoff defining the radius of the sphere
      - the basis vectors of the real space lattice in row
        major format, i.e. a1=cell.row(0), a2=, a3= ...
      -------------------------------------------------------------------
        */

    // Radius of the sphere in reciprocal space defining the
    // maximum spatial resolution
    double kcut = 2.0*PI/(this->smearing);
    // initialization
    math::Kvectors k_vectors(kcut, cell);
    // number of k-vectors
    size_t n_k = k_vectors.get_numvectors();
    // k-vector matrix: the i-th vector is the i-th row of this matrix
    Matrix_t k_vec = k_vectors.get_kvectors();
    // standard norms of the vectors
    Vector_t k_val = k_vectors.get_norms();
    // k_vectors.print_analysis();

    /*
      -------------------------------------------------------------------
      PRECOMPUTE LODE PARAMETERS FROM k-VECTORS:

      Start computing the quantities required by SOAP and/or
      LODE that depend on the wave vectors (and therefore the
      system cell) but NOT on the atomic positions.
      More specifically, these are:
      - the Fourier transformed atomic density:
        G_k = FT[g](\vec{k})
        note that in order to get the LODE representation,
        this is multiplied by the Fourier transformed potential
        as well, so G_k = FT[g]*FT[V]
      - the spherical harmonics Y_lm evaluated at the k-vectors
        Y_lm(\hat{k})
      - the projection of the plane wave basis onto the
        primitive (normalized) radial basis
        I_nl(k) = (R_n(r), j_l(kr))
        ----------------------------------------------------------------
    */

    // Initialize arrays
    // Fourier components of potential / density
    auto G_k{math::Vector_t(n_k)};
    // Spherical harmonics of shape (N_k, (l_max+1)^2)
    auto Y_lm{math::Matrix_t(n_k, n_col)};
    // Radial integrals of shape (N_k, n_max * (l_max + 1)=nl_size) 
    size_t nl_size{this->max_radial*(this->max_angular+1)};
    auto I_nl{math::Matrix_t(n_k, nl_size)};

    // Precompute the three above quantities for each k-vector
    for (size_t ik{0}; ik < n_k; ++ik) {
      // Fourier charge at |k|=k_val (multiply by 1/k^2 for LODE) 
      G_k(ik) = std::exp(-0.5 * pow(this->smearing*k_val(ik) , 2) );

      // Spherical harmonics at hat{k}=vec{k}/k
      if (k_val(ik)!=0.0){
        Eigen::Vector3d k_dir{k_vec.row(ik)/k_val(ik)};
        this->spherical_harmonics.calc(k_dir);
        auto && harmonics{spherical_harmonics.get_harmonics()};
        for (size_t lm{0}; lm < n_col; ++lm) {
          Y_lm(ik, lm) = harmonics(lm);
        }
      } else {
        // For k=0, only Y_00=1/sqrt(4pi) is defined. We divide this
        // by 2 since the k=0 term is not counted twice later on.
        Y_lm(ik, 0) = 0.5 / std::sqrt(4 * PI);
      }

      // Radial integral between primitive radial function and
      // spherical Bessel function (R_n(r), j_l(kr)) up to global
      // factor of sqrt(pi) (irrelevant due to final normalization)
      auto && radint = radial_integral->compute_radial_integral(k_val(ik));
      for (size_t nl{0}; nl < nl_size; ++nl) {
        I_nl(ik, nl) = radint(nl);
      }
    }  // end of loop over k vectors

    /* Test radial scattering function @memo delete after testing phase
    size_t ntest{501};
    double kstep{kcut/(ntest-1)};
    for (size_t ik{0}; ik<ntest; ik++) {
      double kval{kstep * ik};
      std::cout << "k="<<kval<<"\n";
      auto && radint = radial_integral->compute_radial_integral(kval);
      for (size_t nl{0}; nl < nl_size; ++nl){  
        double Inl_test{radint(nl)};
        size_t nrad{nl / (this->max_angular+1)};
        size_t lang{nl % (this->max_angular+1)};
        std::cout << "I_{"<<nrad<<lang<<"} = " << Inl_test << "\n";
      }
      std::cout << "\n";
    } */

    /*
      -------------------------------------------------------------------
      -------------------------------------------------------------------
      ----------START STEPS DEPENDENT ON ATOMIC POSITIONS----------------
      -------------------------------------------------------------------
      -------------------------------------------------------------------
    */

    // Define the basic quantities coming from the atomic positions
    auto tcoords{manager_root->get_positions()};
    Matrix_t coords{tcoords.transpose()};
    size_t natoms = coords.rows();

    // Precompute trigonometric expressions sin(vec{k} * vec{r}), cos(vec{k} * vec{r})
    // for all k-vectors and all atomic positions 
    // Initialization of matrices in which to store results
    auto cos_ki{math::Matrix_t(n_k,natoms)}; 
    auto sin_ki{math::Matrix_t(n_k,natoms)}; 
    for (size_t ik{0}; ik < n_k; ++ik) {
      for (size_t iat{0}; iat < natoms; ++iat) {
        double trigarg{coords.row(iat).dot(k_vec.row(ik))};
        cos_ki(ik,iat) = std::cos(trigarg); 
        sin_ki(ik,iat) = std::sin(trigarg); 
      } 
    }

    /*
    -----------------------------------------------------------------
    COMPLETE EVALUATION OF SPHERICAL EXPANSION:

    All the preparations have been done: this final part combines all
    of the previous results to compute the quantities <anlm|V_i>, the
    spherical expansion coefficients. These will then be used in
    calculate_spherical_invariants to produce quantities invariant
    under rotations.

    This is an O(N^2) implementation in which the two (discrete)
    Fourier transforms are evaluated directly following the
    definition. In the future, we plan to implement an O(NlogN)
    version of the code.
    -----------------------------------------------------------------
    */
    // Global prefactor in expansion coefficients from Fourier transform
    double global_factor = 16.0 * pow(PI,2) / volume;
    // Start the accumulation: loop over all center atoms
    for (auto center : manager) {
      // Preparations:
      auto atom_i_tag{center.get_atom_tag()};
      size_t iat{manager->get_atom_index(atom_i_tag)}; 
      // c^{i}
      auto & coefficients_center = expansions_coefficients[center];
      Key_t center_type{center.get_atom_type()};

      // loop over k-vectors
      for (size_t ik{0}; ik < n_k; ++ik) {
        // Start the accumulation with the central atom contribution
        size_t nl_idx{0};
        for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          size_t l_block_idx{0};
          for (size_t ang_l{0}; ang_l < this->max_angular + 1; ++ang_l) {
            // l odd contributions vanish by k-symmetry for central atom
            size_t size_m{2*ang_l+1};
            if (ang_l%2==0) {
	          for (size_t mval{0}; mval < size_m; ++mval) {
                coefficients_center[center_type](radial_n,l_block_idx+mval) += 
                I_nl(ik,nl_idx) * Y_lm(ik,l_block_idx+mval) * G_k(ik) 
                * global_factor; 
              }
            }
            l_block_idx += size_m;
            nl_idx += 1;
          } // end of loop over l=0,1,...,lmax
        } // end of loop over n=0,1,...,nmax-1

        // Initialize some variables that will be updated in each loop
        double fourier_real{};
        double fourier_imag{};
        double phase_factor{};
        double nlmk_factor{};

        // Start adding up terms for all "neighbor" atoms
        for (auto neigh : manager) {
          // Initialize index and type of new atom and only execute code
          // if different from center atom
          // (the periodic images move with the center, so their
          // contribution to the center gradient is zero as well)
          auto atom_j_tag{neigh.get_atom_tag()};
          size_t jat{manager->get_atom_index(atom_j_tag)};
          if (jat == iat) {
            continue;
          }
          Key_t neigh_type{neigh.get_atom_type()};
          
          // Get real and imaginary parts of the Fourier factor
          // exp(-k*r_ij) that will help us define the phase factor   
          fourier_real = cos_ki(ik,jat)*cos_ki(ik,iat)
            + sin_ki(ik,jat)*sin_ki(ik,iat);
          fourier_imag = sin_ki(ik,jat)*cos_ki(ik,iat)
            - cos_ki(ik,jat)*sin_ki(ik,iat);
 
          // index running over all pairs (n,l), so (0,0)=0, (0,1)=1, etc.
          size_t nl_idx{0};
          for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
            // index running over (l,m) pairs updated in steps of #m=2l+1
            size_t l_block_idx{0};
            for (size_t ang_l{0}; ang_l < this->max_angular+1; ++ang_l) {
              // Get real or imag. part of e^{ikr_ij} based on parity of l
              if (ang_l%2==0) {
                phase_factor = pow(-1.0,ang_l/2) * fourier_real;
              } else {
                phase_factor = -pow(-1.0,(ang_l+1)/2) * fourier_imag;
              }
              size_t size_m = 2*ang_l+1;
              for (size_t mval{0}; mval < size_m; ++mval) {
                size_t lm_idx{l_block_idx + mval};
                // All factors depending on the parameters (n,l,m) and vec{k} 
                nlmk_factor = I_nl(ik,nl_idx) * Y_lm(ik,l_block_idx+mval)
                            * G_k(ik) * global_factor;
                // Add new contribution to expansion coefficient
                coefficients_center[neigh_type](radial_n,l_block_idx+mval) += 
                      phase_factor * nlmk_factor;
          
                if (compute_gradients) {
                  auto & coefficients_center_gradient = 
                     expansions_coefficients_gradient[center.get_atom_ii()];
                  auto && grad_center_by_type{
                     coefficients_center_gradient[neigh_type]};
                  // Get real or imag. part of e^{ikr_ij} based on parity of l
                  if (ang_l%2==0) {
                    phase_factor = pow(-1.0,ang_l/2) * fourier_imag;
                  } else {
                    phase_factor = -pow(-1.0,(ang_l+1)/2) * fourier_real;
                  }
                  // Update x,y and z components of gradients
                  for (size_t car_idx{0}; car_idx < ThreeD; car_idx++) {
                  grad_center_by_type(car_idx*max_radial+radial_n, lm_idx) +=
                      phase_factor * nlmk_factor * k_vec(ik,car_idx);
                  } // for (cartesian_idx)
                } // if (compute_gradients)
              } // for (mval)
              l_block_idx += size_m;
              nl_idx += 1;
            } // for (ang_l)
          } // for (radial_n)         
        } // for (neigh : center)
      } // for (kvectors) 

      // Normalize and orthogonalize the radial coefficients
      radial_integral->finalize_coefficients(coefficients_center);
      
      // Write code in which to store the coefficients 
      std::ofstream myfile;
      myfile.open ("expansioncoefficients_kspace.txt", std::ios::app);
      for (auto neigh : manager) { 
        auto atom_j_tag{neigh.get_atom_tag()};
        size_t jat{manager->get_atom_index(atom_j_tag)};
        Key_t neigh_type{neigh.get_atom_type()};
        auto && coeff = coefficients_center[neigh_type];
        myfile << "Indices " << iat << " and " << jat << "\n";
        myfile << coeff.rows() << " x " <<
               coeff.cols() << " coeff matrix = \n" <<
               coeff << "\n";
      }
      myfile.close();

      // if (compute_gradients) {
      //  radial_integral->template finalize_coefficients_der<ThreeD>(
      //      expansions_coefficients_gradient, center);
      //}
    } // for (center : manager)
  } // compute()

  // STRUCTURE MANAGER STUFF BELOW
  template <class StructureManager>
  void
  CalculatorKspaceSphericalExpansion::initialize_expansion_environment_wise(
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
  void
  CalculatorKspaceSphericalExpansion::initialize_expansion_with_global_species(
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
    static rascal::CalculatorKspaceSphericalExpansion
    from_json(const json & j) {
      return rascal::CalculatorKspaceSphericalExpansion{j};
    }

    static void to_json(json & j,
                        const rascal::CalculatorKspaceSphericalExpansion & t) {
      j = t.hypers;
    }
  };
}  // namespace nlohmann

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_KSPACE_HH_
