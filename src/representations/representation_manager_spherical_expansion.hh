/**
 * file   representation_manager_spherical_expansion.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Andrea Grifasi <andrea.grifasi@epfl.ch>
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

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_HH_

#include "representations/representation_manager_base.hh"
#include "representations/cutoff_functions.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include "math/math_utils.hh"
#include "structure_managers/property_block_sparse.hh"

#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace rascal {

  struct AtomicSmearingSpecificationBase {
    //! Constructor
    AtomicSmearingSpecificationBase() = default;
    //! Destructor
    virtual ~AtomicSmearingSpecificationBase() = default;
    //! Copy constructor
    AtomicSmearingSpecificationBase(
        const AtomicSmearingSpecificationBase & other) = delete;
    //! Move constructor
    AtomicSmearingSpecificationBase(AtomicSmearingSpecificationBase && other) =
        default;
    //! Copy assignment operator
    AtomicSmearingSpecificationBase &
    operator=(const AtomicSmearingSpecificationBase & other) = delete;
    //! Move assignment operator
    AtomicSmearingSpecificationBase &
    operator=(AtomicSmearingSpecificationBase && other) = default;

    using Hypers_t = RepresentationManagerBase::Hypers_t;
  };

  namespace internal {

    enum class GaussianSigmaType { Constant, PerSpecies, Radial, End_ };

    /**
     * Specification to hold the parameter for the atomic smearing function,
     * currently only Gaussians are supported.
     *
     * This is `sigma' in the definition `f(r) = A exp(r / (2 sigma^2))'.
     * The width may depend both on the atomic species of the neighbour as well
     * as the distance.
     *
     * Note that this function is template-specialized by Gaussian sigma type
     * (constant, per-species, or radially dependent).
     *
     * @param pair Atom pair defining the neighbour, as e.g. returned by
     *             iteration over neighbours of a centre
     *
     * @throw logic_error if the requested sigma type has not been implemented
     *
     */
    template <GaussianSigmaType SigmaType>
    struct AtomicSmearingSpecification {};

    template <>
    struct AtomicSmearingSpecification<GaussianSigmaType::Constant>
        : AtomicSmearingSpecificationBase {
      using Hypers_t = typename AtomicSmearingSpecificationBase::Hypers_t;
      explicit AtomicSmearingSpecification(const Hypers_t & hypers) {
        this->constant_gaussian_sigma =
            hypers.at("gaussian_sigma").at("value").get<double>();
      }
      template <size_t Order, size_t Layer>
      double get_gaussian_sigma(ClusterRefKey<Order, Layer> & /* pair */) {
        return this->constant_gaussian_sigma;
      }
      double constant_gaussian_sigma{0.};
    };

    /** Per-species template specialization of the above */

    template <>
    struct AtomicSmearingSpecification<GaussianSigmaType::PerSpecies>
        : AtomicSmearingSpecificationBase {
      using Hypers_t = typename AtomicSmearingSpecificationBase::Hypers_t;
      explicit AtomicSmearingSpecification(const Hypers_t & /* hypers */) {}
      template <size_t Order, size_t Layer>
      double get_gaussian_sigma(ClusterRefKey<Order, Layer> & /* pair */) {
        throw std::logic_error("Requested a sigma type that has not yet "
                               "been implemented");
        return -1;
      }
    };

    /** Radially-dependent template specialization of the above */
    template <>
    struct AtomicSmearingSpecification<GaussianSigmaType::Radial>
        : AtomicSmearingSpecificationBase {
      using Hypers_t = typename AtomicSmearingSpecificationBase::Hypers_t;
      explicit AtomicSmearingSpecification(const Hypers_t & /* hypers */) {}
      template <size_t Order, size_t Layer>
      double get_gaussian_sigma(ClusterRefKey<Order, Layer> & /* pair */) {
        throw std::logic_error("Requested a sigma type that has not yet "
                               "been implemented");
        return -1;
      }
    };

  }  // namespace internal

  template <internal::GaussianSigmaType Type, class Hypers>
  decltype(auto) make_gaussian_density(const Hypers & sigma_hypers) {
    return std::static_pointer_cast<AtomicSmearingSpecificationBase>(
        std::make_shared<internal::AtomicSmearingSpecification<Type>>(
            sigma_hypers));
  }

  template <internal::GaussianSigmaType Type>
  decltype(auto) downcast_gaussian_density(
      const std::shared_ptr<AtomicSmearingSpecificationBase> &
          gaussian_density) {
    return std::static_pointer_cast<
        internal::AtomicSmearingSpecification<Type>>(gaussian_density);
  }

  struct RadialContributionBase {
    //! Constructor
    RadialContributionBase() = default;
    //! Destructor
    virtual ~RadialContributionBase() = default;
    //! Copy constructor
    RadialContributionBase(const RadialContributionBase & other) = delete;
    //! Move constructor
    RadialContributionBase(RadialContributionBase && other) = default;
    //! Copy assignment operator
    RadialContributionBase &
    operator=(const RadialContributionBase & other) = delete;
    //! Move assignment operator
    RadialContributionBase &
    operator=(RadialContributionBase && other) = default;

    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using Matrix_t = Eigen::MatrixXd;
    using Vector_t = Eigen::VectorXd;
    using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;
    using Vector_Ref = typename Eigen::Ref<const Vector_t>;

    //! Pure Virtual Function to set hyperparameters of the cutoff function
    virtual void set_hyperparameters(const Hypers_t &) = 0;

    virtual void precompute() = 0;
    virtual Vector_Ref compute_center_contribution(const double & fac_a) = 0;
    virtual Matrix_Ref compute_neighbour_contribution(const double & distance,
                                                      const double & fac_a) = 0;
  };

  struct RadialContributionGTO : RadialContributionBase {
    //! Constructor
    explicit RadialContributionGTO(const Hypers_t & hypers) {
      this->set_hyperparameters(hypers);
      this->precompute();
    }
    //! Destructor
    virtual ~RadialContributionGTO() = default;
    //! Copy constructor
    RadialContributionGTO(const RadialContributionGTO & other) = delete;
    //! Move constructor
    RadialContributionGTO(RadialContributionGTO && other) = default;
    //! Copy assignment operator
    RadialContributionGTO &
    operator=(const RadialContributionGTO & other) = delete;
    //! Move assignment operator
    RadialContributionGTO & operator=(RadialContributionGTO && other) = default;

    using Parent = RadialContributionBase;
    using Hypers_t = typename Parent::Hypers_t;
    using Matrix_t = typename Parent::Matrix_t;
    using Vector_t = typename Parent::Vector_t;
    using Matrix_Ref = typename Parent::Matrix_Ref;
    using Vector_Ref = typename Parent::Vector_Ref;

    //! Pure Virtual Function to set hyperparameters of the cutoff function
    void set_hyperparameters(const Hypers_t & hypers) {
      this->hypers = hypers;

      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");

      this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
      this->fac_b.resize(this->max_radial, 1);
      this->radial_norm_factors.resize(this->max_radial, 1);
      this->radial_nl_factors.resize(this->max_radial, this->max_angular + 1);
      this->radial_sigmas.resize(this->max_radial, 1);
      this->radial_integral_neighbour.resize(this->max_radial,
                                             this->max_angular + 1);
      this->radial_integral_center.resize(this->max_radial);

      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
    }

    void precompute() {
      this->precompute_radial_sigmas();
      this->precompute_radial_overlap();
    }

    Vector_Ref compute_center_contribution(const double & fac_a) {
      using math::PI;
      using std::pow;
      using std::sqrt;
      // Contribution from the central atom
      // Y_l^m(θ, φ) =  √((2l+1)/(4*π))) \delta_{m0} and
      //  \delta_{l0} (spherical symmetry) and
      // 1F1(a,b,0) = 1
      for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
        this->radial_integral_center(radial_n) =
            this->radial_norm_factors(radial_n) *
            this->radial_nl_factors(radial_n, 0) * 0.25 *
            pow(fac_a + this->fac_b[radial_n], -0.5 * (3.0 + radial_n));
      }

      this->radial_integral_center =
          this->radial_ortho_matrix * this->radial_integral_center;
      return Vector_Ref(this->radial_integral_center);
    }

    Matrix_Ref compute_neighbour_contribution(const double & distance,
                                              const double & fac_a) {
      double exp_factor = std::exp(-fac_a * pow(distance, 2));
      for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
        double radial_sigma_factors{pow(fac_a, 2) /
                                    (fac_a + this->fac_b[radial_n])};
        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             angular_l++) {
          this->radial_integral_neighbour(radial_n, angular_l) =
              exp_factor * this->radial_nl_factors(radial_n, angular_l) * 0.25 *
              pow(fac_a + this->fac_b[radial_n],
                  -0.5 * (3.0 + angular_l + radial_n)) *
              pow(distance * fac_a, angular_l) *
              math::hyp1f1(0.5 * (3.0 + angular_l + radial_n), 1.5 + angular_l,
                           pow(distance, 2) * radial_sigma_factors);
        }
        this->radial_integral_neighbour.row(radial_n) *=
            this->radial_norm_factors(radial_n);
      }
      this->radial_integral_neighbour =
          this->radial_ortho_matrix * this->radial_integral_neighbour;
      return Matrix_Ref(this->radial_integral_neighbour);
    }

    /** Compute common prefactors for the radial Gaussian basis functions */
    void precompute_radial_sigmas() {
      using std::pow;

      for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
        this->radial_sigmas[radial_n] =
            std::max(std::sqrt(static_cast<double>(radial_n)), 1.0) *
            this->interaction_cutoff / static_cast<double>(this->max_radial);
        this->fac_b[radial_n] = 0.5 * pow(this->radial_sigmas[radial_n], -2);
      }

      // Precompute common prefactors
      for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
        this->radial_norm_factors(radial_n) = std::sqrt(
            2.0 / (std::tgamma(1.5 + radial_n) *
                   pow(this->radial_sigmas[radial_n], 3.0 + 2.0 * radial_n)));
        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             ++angular_l) {
          this->radial_nl_factors(radial_n, angular_l) =
              std::tgamma(0.5 * (3.0 + angular_l + radial_n)) /
              std::tgamma(1.5 + angular_l);
        }
      }
    }

    /**
     * Compute the radial overlap matrix for later orthogonalization.
     *
     * @throw runtime_error if the overlap matrix cannot be diagonalized
     */
    void precompute_radial_overlap() {
      using std::pow;
      using std::sqrt;
      using std::tgamma;

      // TODO(max-veit) see if we can replace the gammas with their natural
      // logs,
      // since it'll overflow for relatively small n (n1 + n2 >~ 300)
      Eigen::MatrixXd overlap(this->max_radial, this->max_radial);
      for (size_t radial_n1{0}; radial_n1 < this->max_radial; radial_n1++) {
        for (size_t radial_n2{0}; radial_n2 < this->max_radial; radial_n2++) {
          overlap(radial_n1, radial_n2) =
              pow(0.5 / pow(this->radial_sigmas[radial_n1], 2) +
                      0.5 / pow(this->radial_sigmas[radial_n2], 2),
                  -0.5 * (3.0 + radial_n1 + radial_n2)) /
              (pow(this->radial_sigmas[radial_n1], radial_n1) *
               pow(this->radial_sigmas[radial_n2], radial_n2)) *
              tgamma(0.5 * (3.0 + radial_n1 + radial_n2)) /
              (pow(this->radial_sigmas[radial_n1] *
                       this->radial_sigmas[radial_n2],
                   1.5) *
               sqrt(tgamma(1.5 + radial_n1) * tgamma(1.5 + radial_n2)));
        }
      }

      // Compute the inverse square root of the overlap matrix
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(overlap);
      if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Warning: Could not diagonalize "
                                 "radial overlap matrix.");
      }
      Eigen::MatrixXd eigenvalues = eigensolver.eigenvalues();
      Eigen::ArrayXd eigs_invsqrt = eigenvalues.array().sqrt().inverse();
      Eigen::MatrixXd unitary = eigensolver.eigenvectors();
      this->radial_ortho_matrix =
          unitary * eigs_invsqrt.matrix().asDiagonal() * unitary.adjoint();
    }

    Matrix_t radial_integral_neighbour{};
    Vector_t radial_integral_center{};

    Hypers_t hypers{};
    double interaction_cutoff{};
    size_t max_radial{};
    size_t max_angular{};
    size_t n_species{};

    Vector_t radial_sigmas{};
    // b = 1 / (2*\sigma_n^2)
    Vector_t fac_b{};
    Vector_t radial_norm_factors{};
    Vector_t radial_sigma_factors{};
    Matrix_t radial_nl_factors{};
    Matrix_t radial_ortho_matrix{};
  };

  /**
   * Handles the expansion of an environment in a spherical and radial basis.
   *
   * The local environment of each atom is represented by Gaussians of a
   * certain width (user-defined; can be constant, species-dependent, or
   * radially dependent).  This density field is expanded in an angular basis
   * of spherical harmonics (à la SOAP) and a radial basis of either Gaussians
   * (again, as in SOAP) or one of the more recent bases currently under
   * development.
   */

  template <class StructureManager>
  class RepresentationManagerSphericalExpansion
      : public RepresentationManagerBase {
   public:
    using Parent = RepresentationManagerBase;
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = typename Parent::Hypers_t;
    using ReferenceHypers_t = Parent::ReferenceHypers_t;
    using Key_t = std::vector<int>;
    using SparseProperty_t = BlockSparseProperty<double, 1, 0>;
    using Dense_t = typename SparseProperty_t::Dense_t;
    using Data_t = typename SparseProperty_t::Data_t;
    /**
     * Set the hyperparameters of this descriptor from a json object.
     *
     * @param hypers structure (usually parsed from json) containing the
     *               options and hyperparameters
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the structure
     */
    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::CutoffFunctionType;
      using internal::GaussianSigmaType;

      this->hypers = hypers;

      this->expansions_coefficients.clear();

      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      if (hypers.find("n_species") != hypers.end()) {
        this->n_species = hypers.at("n_species");
      } else {
        this->n_species = 1;  // default: no species distinction
      }

      auto radial_contribution_hypers =
          hypers.at("radial_contribution").get<json>();
      auto radial_contribution_type =
          radial_contribution_hypers.at("type").get<std::string>();
      if (radial_contribution_type.compare("GTO") == 0) {
        this->radial_integral =
            std::static_pointer_cast<RadialContributionBase>(
                std::make_shared<RadialContributionGTO>(hypers));
      } else {
        throw std::logic_error(
            "Requested Radial contribution type \'" + radial_contribution_type +
            "\' has not been implemented.  Must be one of" + ": \'GTO\'.");
      }

      auto density_hypers = hypers.at("gaussian_density").get<json>();
      auto density_type = density_hypers.at("type").get<std::string>();
      if (density_type.compare("Constant") == 0) {
        this->gaussian_density_type = GaussianSigmaType::Constant;
        this->gaussian_density =
            make_gaussian_density<GaussianSigmaType::Constant>(density_hypers);
      } else {
        throw std::logic_error(
            "Requested Gaussian sigma type \'" + density_type +
            "\' has not been implemented.  Must be one of" + ": \'Constant\'.");
      }

      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      auto fc_type = fc_hypers.at("type").get<std::string>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
      this->cutoff_smooth_width = fc_hypers.at("smooth_width").at("value");
      if (fc_type.compare("Cosine") == 0) {
        this->cutoff_function_type = CutoffFunctionType::Cosine;
        this->cutoff_function =
            make_cutoff_function<CutoffFunctionType::Cosine>(fc_hypers);
      } else {
        throw std::logic_error("Requested cutoff function type \'" + fc_type +
                               "\' has not been implemented.  Must be one of" +
                               ": \'Cosine\'.");
      }
    }

    /**
     * Construct a new RepresentationManager using a hyperparameters container
     *
     * @param hypers container (usually parsed from json) for the options and
     *               hyperparameters
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the container
     */
    RepresentationManagerSphericalExpansion(ManagerPtr_t sm,
                                            const Hypers_t & hyper)
        : expansions_coefficients{*sm}, structure_manager{std::move(sm)} {
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    RepresentationManagerSphericalExpansion(
        const RepresentationManagerSphericalExpansion & other) = delete;

    //! Move constructor
    RepresentationManagerSphericalExpansion(
        RepresentationManagerSphericalExpansion && other) = default;

    //! Destructor
    virtual ~RepresentationManagerSphericalExpansion() = default;

    //! Copy assignment operator
    RepresentationManagerSphericalExpansion &
    operator=(const RepresentationManagerSphericalExpansion & other) = delete;

    //! Move assignment operator
    RepresentationManagerSphericalExpansion &
    operator=(RepresentationManagerSphericalExpansion && other) = default;

    //! compute representation
    void compute();

    template <internal::GaussianSigmaType SigmaType,
              internal::CutoffFunctionType FcType>
    void compute_by_gaussian_density(
        const std::shared_ptr<internal::CutoffFunction<FcType>> & fc,
        const std::shared_ptr<
            internal::AtomicSmearingSpecification<SigmaType>> & sigma);

    std::vector<Precision_t> & get_representation_raw_data() {
      return this->dummy;
    }

    Data_t & get_representation_sparse_raw_data() {
      return this->expansions_coefficients.get_raw_data();
    }

    size_t get_feature_size() {
      return this->expansions_coefficients.get_nb_comp();
    }

    size_t get_center_size() {
      return this->expansions_coefficients.get_nb_item();
    }

    auto get_representation_full() {
      return this->expansions_coefficients.get_dense_rep();
    }

    SparseProperty_t expansions_coefficients;

   protected:
   private:
    double interaction_cutoff{};
    double cutoff_smooth_width{};
    size_t max_radial{};
    size_t max_angular{};
    size_t n_species{};

    std::vector<Precision_t> dummy{};

    ManagerPtr_t structure_manager;

    internal::GaussianSigmaType gaussian_density_type{};
    internal::CutoffFunctionType cutoff_function_type{};

    std::shared_ptr<RadialContributionBase> radial_integral{};
    std::shared_ptr<CutoffFunctionBase> cutoff_function{};
    std::shared_ptr<AtomicSmearingSpecificationBase> gaussian_density{};

    Hypers_t hypers{};
  };

  template <class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::compute() {
    using internal::CutoffFunctionType;
    using internal::GaussianSigmaType;

    internal::CombineEnums<GaussianSigmaType, CutoffFunctionType> combineEnums;
    switch (
        combineEnums(this->gaussian_density_type, this->cutoff_function_type)) {
    case combineEnums(GaussianSigmaType::Constant,
                      CutoffFunctionType::Cosine): {
      auto fc{downcast_cutoff_function<CutoffFunctionType::Cosine>(
          this->cutoff_function)};
      auto sigma{downcast_gaussian_density<GaussianSigmaType::Constant>(
          this->gaussian_density)};
      this->compute_by_gaussian_density(fc, sigma);
      break;
    }
    case combineEnums(GaussianSigmaType::PerSpecies,
                      CutoffFunctionType::Cosine): {
      auto fc{downcast_cutoff_function<CutoffFunctionType::Cosine>(
          this->cutoff_function)};
      auto sigma{downcast_gaussian_density<GaussianSigmaType::PerSpecies>(
          this->gaussian_density)};
      this->compute_by_gaussian_density(fc, sigma);
      break;
    }
    case combineEnums(GaussianSigmaType::Radial, CutoffFunctionType::Cosine): {
      auto fc{downcast_cutoff_function<CutoffFunctionType::Cosine>(
          this->cutoff_function)};
      auto sigma{downcast_gaussian_density<GaussianSigmaType::Radial>(
          this->gaussian_density)};
      this->compute_by_gaussian_density(fc, sigma);
      break;
    }

    default:
      throw std::logic_error("The combination of parameter is not handdled.");
      break;
    }
  }

  /**
   * Compute the spherical expansion
   * TODO(felix,max) use the parity of the spherical harmonics to use half
   * neighbourlist, i.e. C^{ij}_{nlm} = (-1)^l C^{ji}_{nlm}.
   */
  template <class Mngr>
  template <internal::GaussianSigmaType SigmaType,
            internal::CutoffFunctionType FcType>
  void
  RepresentationManagerSphericalExpansion<Mngr>::compute_by_gaussian_density(
      const std::shared_ptr<internal::CutoffFunction<FcType>> & cutoff_function,
      const std::shared_ptr<internal::AtomicSmearingSpecification<SigmaType>> &
          gaussian_density) {
    using math::PI;
    using std::pow;

    auto n_row{this->max_radial};
    auto n_col{pow(this->max_angular + 1, 2)};
    this->expansions_coefficients.clear();
    this->expansions_coefficients.set_shape(n_row, n_col);
    this->expansions_coefficients.resize();

    for (auto center : this->structure_manager) {
      auto & coefficients_center = this->expansions_coefficients[center];

      // initialize the expansion coefficients to 0
      Key_t center_type{center.get_atom_type()};
      coefficients_center[center_type] = Dense_t::Zero(n_row, n_col);
      for (auto neigh : center) {
        Key_t neigh_type{neigh.get_atom_type()};
        // avoid initializing again same chemical channel
        if (coefficients_center.count(neigh_type) == 0) {
          coefficients_center[neigh_type] = Dense_t::Zero(n_row, n_col);
        }
      }

      // Eigen::MatrixXd radial_integral =
      //     Eigen::MatrixXd::Zero(this->max_radial, this->max_angular + 1);

      // Start the accumulator with the central atom
      // All terms where l =/= 0 cancel

      // a = 1 / (2*\sigma^2)
      double fac_a{0.5 * pow(gaussian_density->get_gaussian_sigma(center), -2)};

      // // Contribution from the central atom
      // // Y_l^m(θ, φ) =  √((2l+1)/(4*π))) \delta_{m0} and
      // //  \delta_{l0} (spherical symmetry) and
      // // 1F1(a,b,0) = 1
      // for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
      //   radial_integral(radial_n, 0) =
      //       radial_norm_factors(radial_n) * radial_nl_factors(radial_n, 0) *
      //       0.25 * pow(fac_a + this->fac_b[radial_n], -0.5 * (3.0 +
      //       radial_n));
      // }
      // coefficients_center[center_type].col(0) +=
      //     this->radial_ortho_matrix * radial_integral.col(0) / sqrt(4.0 *
      //     PI);
      coefficients_center[center_type].col(0) +=
          this->radial_integral->compute_center_contribution(fac_a) /
          sqrt(4.0 * PI);

      for (auto neigh : center) {
        auto dist{this->structure_manager->get_distance(neigh)};
        auto direction{this->structure_manager->get_direction_vector(neigh)};
        fac_a = 0.5 * pow(gaussian_density->get_gaussian_sigma(neigh), -2);

        Key_t neigh_type{neigh.get_atom_type()};

        // Note: the copy _should_ be optimized out (RVO)
        Eigen::MatrixXd harmonics =
            math::compute_spherical_harmonics(direction, this->max_angular);

        // for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
        //   // TODO(max-veit) this is all specific to the Gaussian radial basis
        //   // (doing just the angular integration would instead spit out
        //   // spherical Bessel functions below)

        //   double radial_sigma_factors{pow(fac_a, 2) /
        //                               (fac_a + this->fac_b[radial_n])};
        //   for (size_t angular_l{0}; angular_l < this->max_angular + 1;
        //        angular_l++) {
        //     radial_integral(radial_n, angular_l) =
        //         exp_factor * radial_nl_factors(radial_n, angular_l) * 0.25 *
        //         pow(fac_a + this->fac_b[radial_n],
        //             -0.5 * (3.0 + angular_l + radial_n)) *
        //         pow(dist * fac_a, angular_l) *
        //         math::hyp1f1(0.5 * (3.0 + angular_l + radial_n),
        //                      1.5 + angular_l,
        //                      pow(dist, 2) * radial_sigma_factors);
        //   }
        //   radial_integral.row(radial_n) *= radial_norm_factors(radial_n);
        // }
        // radial_integral = this->radial_ortho_matrix * radial_integral;

        auto neighbour_contribution =
            this->radial_integral->compute_neighbour_contribution(dist, fac_a);
        harmonics *= cutoff_function->f_c(dist);
        for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
          size_t lm_collective_idx{0};
          for (size_t angular_l{0}; angular_l < this->max_angular + 1;
               angular_l++) {
            for (size_t m_array_idx{0}; m_array_idx < 2 * angular_l + 1;
                 m_array_idx++) {
              coefficients_center[neigh_type](radial_n, lm_collective_idx) +=
                  neighbour_contribution(radial_n, angular_l) *
                  harmonics(angular_l, m_array_idx);
              ++lm_collective_idx;
            }
          }
        }
      }  // for (neigh : center)
    }    // for (center : structure_manager)
  }      // compute()

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_HH_
