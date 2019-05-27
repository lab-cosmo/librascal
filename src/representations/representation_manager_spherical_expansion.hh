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
#include "math/hyp1f1.hh"
#include "structure_managers/property_block_sparse.hh"

#include <algorithm>
#include <cmath>
#include <memory>
#include <exception>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unordered_set>

namespace rascal {

  namespace internal {

    /**
     * List of possible Radial basis that can be used by the spherical
     * expansion.
     */
    enum class RadialBasisType { GTO, End_ };

    /**
     * List of possible atomic smearing for the definition of the atomic
     * density. If not specified, the gaussian type smearing is implided.
     */
    enum class AtomicSmearingType { Constant, PerSpecies, Radial, End_ };

    /**
     * Base class for the specification of the atomic smearing.
     */
    struct AtomicSmearingSpecificationBase {
      //! Constructor
      AtomicSmearingSpecificationBase() = default;
      //! Destructor
      virtual ~AtomicSmearingSpecificationBase() = default;
      //! Copy constructor
      AtomicSmearingSpecificationBase(
          const AtomicSmearingSpecificationBase & other) = delete;
      //! Move constructor
      AtomicSmearingSpecificationBase(
          AtomicSmearingSpecificationBase && other) = default;
      //! Copy assignment operator
      AtomicSmearingSpecificationBase &
      operator=(const AtomicSmearingSpecificationBase & other) = delete;
      //! Move assignment operator
      AtomicSmearingSpecificationBase &
      operator=(AtomicSmearingSpecificationBase && other) = default;

      using Hypers_t = RepresentationManagerBase::Hypers_t;
    };

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
    template <AtomicSmearingType SigmaType>
    struct AtomicSmearingSpecification {};

    template <>
    struct AtomicSmearingSpecification<AtomicSmearingType::Constant>
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
    struct AtomicSmearingSpecification<AtomicSmearingType::PerSpecies>
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
    struct AtomicSmearingSpecification<AtomicSmearingType::Radial>
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

    //! Utility to make shared pointer and cast to base class
    template <AtomicSmearingType Type, class Hypers>
    decltype(auto) make_atomic_smearing(const Hypers & sigma_hypers) {
      return std::static_pointer_cast<AtomicSmearingSpecificationBase>(
          std::make_shared<AtomicSmearingSpecification<Type>>(sigma_hypers));
    }

    //! Utility to cast base to child class
    template <AtomicSmearingType Type>
    decltype(auto) downcast_atomic_smearing(
        const std::shared_ptr<AtomicSmearingSpecificationBase> &
            atomic_smearing) {
      return std::static_pointer_cast<AtomicSmearingSpecification<Type>>(
          atomic_smearing);
    }

    /**
     * Base class to define the radial contribution to the spherical expansion
     */
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
      using Matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;
      using Vector_t = Eigen::VectorXd;
      using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;
      using Vector_Ref = typename Eigen::Ref<const Vector_t>;

      //! Pure Virtual Function to set hyperparameters of the cutoff function
      virtual void set_hyperparameters(const Hypers_t &) = 0;

      virtual void precompute() = 0;
      // Can't make templated virtual member function... But these functions
      // are expected
      // virtual Vector_Ref compute_center_contribution() = 0;
      // virtual Matrix_Ref compute_neighbour_contribution() = 0;
    };

    template <RadialBasisType RBT>
    struct RadialContribution {};

    /**
     * Implementation of the radial contribution for Gaussian Type Orbitals
     * radial basis functions and gaussian smearing of the atom density.
     *
     *          R^{GTO}_{n}(r) = \mathcal{N}_n\ r^{n} \exp[-br^2]
     *
     * \mathcal{N}_n^2 = \frac{2}{\sigma_n^{2n + 3}\Gamma(n + 3/2)}
     * \sigma_n = (r_\text{cut}-\delta r_\text{cut})
     * \max(\sqrt{n},1)/n_\text{max} b=\frac{1}{2\sigma_n} \int_0^\infty
     * R^{GTO}_{n}(r) R^{GTO}_{n\prime}(r) \dd{r}= 2 \left(\frac{1}{2
     * \sigma_{n}^2}+\frac{1}{2 \sigma_{n\prime}^2} \right)^{-\frac{1}{2}
     * (3+n+n\prime)} \Gamma(\frac{3+n+n\prime}{2})
     */
    template <>
    struct RadialContribution<RadialBasisType::GTO> : RadialContributionBase {
      //! Constructor
      explicit RadialContribution(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
        this->precompute();
      }
      //! Destructor
      virtual ~RadialContribution() = default;
      //! Copy constructor
      RadialContribution(const RadialContribution & other) = delete;
      //! Move constructor
      RadialContribution(RadialContribution && other) = default;
      //! Copy assignment operator
      RadialContribution & operator=(const RadialContribution & other) = delete;
      //! Move assignment operator
      RadialContribution & operator=(RadialContribution && other) = default;

      using Parent = RadialContributionBase;
      using Hypers_t = typename Parent::Hypers_t;
      using Matrix_t = typename Parent::Matrix_t;
      using Vector_t = typename Parent::Vector_t;
      using Matrix_Ref = typename Parent::Matrix_Ref;
      using Vector_Ref = typename Parent::Vector_Ref;

      /**
       * Set hyperparameters.
       * @params hypers is expected to be the same as the the input of
       *         the spherical expansion
       */
      void set_hyperparameters(const Hypers_t & hypers) {
        this->hypers = hypers;

        this->max_radial = hypers.at("max_radial");
        this->max_angular = hypers.at("max_angular");

        // init size of the member data
        // both precomputed quantities and actual expansion coefficients
        this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
        this->fac_b.resize(this->max_radial, 1);
        this->radial_norm_factors.resize(this->max_radial, 1);
        this->radial_n_factors.resize(this->max_radial);
        this->radial_sigmas.resize(this->max_radial, 1);
        this->radial_integral_neighbour.resize(this->max_radial,
                                               this->max_angular + 1);
        this->radial_integral_center.resize(this->max_radial);

        // find the cutoff radius of the representation
        auto fc_hypers = hypers.at("cutoff_function").get<json>();
        this->interaction_cutoff =
            fc_hypers.at("cutoff").at("value").get<double>();

        // define the type of smearing to use
        auto smearing_hypers = hypers.at("gaussian_density").get<json>();
        auto smearing_type = smearing_hypers.at("type").get<std::string>();
        if (smearing_type.compare("Constant") == 0) {
          this->atomic_smearing_type = AtomicSmearingType::Constant;
          this->atomic_smearing =
              make_atomic_smearing<AtomicSmearingType::Constant>(
                  smearing_hypers);
        } else {
          throw std::logic_error(
              "Requested Gaussian sigma type \'" + smearing_type +
              "\' has not been implemented.  Must be one of" +
              ": \'Constant\'.");
        }
      }

      void precompute() {
        this->precompute_radial_sigmas();
        this->precompute_radial_overlap();
        this->hyp1f1_calculator.precompute(this->max_radial, this->max_angular);
      }

      //! define the contribution from the central atom to the expansion
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Vector_Ref
      compute_center_contribution(ClusterRefKey<Order, Layer> & center) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};

        // a = 1 / (2*\sigma^2)
        double fac_a{0.5 * pow(smearing->get_gaussian_sigma(center), -2)};

        // Contribution from the central atom
        // Y_l^m(θ, φ) =  √((2l+1)/(4*π))) \delta_{m0} and
        //  \delta_{l0} (spherical symmetry) and
        // 1F1(a,b,0) = 1
        for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
          this->radial_integral_center(radial_n) =
              this->radial_norm_factors(radial_n) *
              this->radial_n_factors(radial_n) *
              sqrt(pow(fac_a + this->fac_b[radial_n], -(3 + radial_n)));
        }

        this->radial_integral_center =
            this->radial_ortho_matrix * this->radial_integral_center;
        return Vector_Ref(this->radial_integral_center);
      }

      //! define the contribution from a neighbour atom to the expansion
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(const double & distance,
                                     ClusterRefKey<Order, Layer> & pair) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};
        // a = 1 / (2*\sigma^2)
        double fac_a{0.5 * pow(smearing->get_gaussian_sigma(pair), -2)};

        // computes (r_{ij}*a)^l incrementally
        Vector_t distance_fac_a_l(this->max_angular + 1);
        distance_fac_a_l(0) = 1.;
        double distance_fac_a{distance * fac_a};
        for (size_t angular_l{1}; angular_l < this->max_angular + 1;
             angular_l++) {
          distance_fac_a_l(angular_l) =
              distance_fac_a_l(angular_l - 1) * distance_fac_a;
        }

        // computes (a+b_n)^{-0.5*(3+l+n)}
        Matrix_t a_b_l_n(this->max_radial, this->max_angular + 1);
        for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
          double a_b_l{1. / sqrt(fac_a + this->fac_b[radial_n])};
          a_b_l_n(radial_n, 0) =
              sqrt(pow(fac_a + this->fac_b[radial_n], -(3 + radial_n)));
          for (size_t angular_l{1}; angular_l < this->max_angular + 1;
               angular_l++) {
            a_b_l_n(radial_n, angular_l) =
                a_b_l_n(radial_n, angular_l - 1) * a_b_l;
          }
        }

        this->hyp1f1_calculator.calc(distance, fac_a, this->fac_b);

        this->radial_integral_neighbour =
            (a_b_l_n.array() * this->hyp1f1_calculator.get_values().array())
                .matrix() * distance_fac_a_l.asDiagonal();
        this->radial_integral_neighbour.transpose() *=
            this->radial_norm_factors.asDiagonal();
        this->radial_integral_neighbour.transpose() *=
            this->radial_ortho_matrix;
        return Matrix_Ref(this->radial_integral_neighbour);
      }

      /** Compute common prefactors for the radial Gaussian basis functions */
      void precompute_radial_sigmas() {
        using math::pow;
        // Precompute common prefactors
        for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          this->radial_sigmas[radial_n] =
              std::max(std::sqrt(static_cast<double>(radial_n)), 1.0) *
              this->interaction_cutoff / static_cast<double>(this->max_radial);
          this->fac_b[radial_n] = 0.5 * pow(this->radial_sigmas[radial_n], -2);
          this->radial_norm_factors(radial_n) =
              0.25 * std::sqrt(2.0 / (std::tgamma(1.5 + radial_n) *
                                      pow(this->radial_sigmas[radial_n],
                                          3 + 2 * radial_n)));
          this->radial_n_factors(radial_n) =
              std::tgamma(0.5 * (3.0 + radial_n)) / std::tgamma(1.5);
        }
      }

      /**
       * Compute the radial overlap matrix for later orthogonalization.
       *
       * @throw runtime_error if the overlap matrix cannot be diagonalized
       */
      void precompute_radial_overlap() {
        using math::pow;
        using std::sqrt;
        using std::tgamma;

        // TODO(max-veit) see if we can replace the gammas with their natural
        // logs,
        // since it'll overflow for relatively small n (n1 + n2 >~ 300)
        Matrix_t overlap(this->max_radial, this->max_radial);
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
        Eigen::SelfAdjointEigenSolver<Matrix_t> eigensolver(overlap);
        if (eigensolver.info() != Eigen::Success) {
          throw std::runtime_error("Warning: Could not diagonalize "
                                   "radial overlap matrix.");
        }
        Matrix_t eigenvalues = eigensolver.eigenvalues();
        Eigen::ArrayXd eigs_invsqrt = eigenvalues.array().sqrt().inverse();
        Matrix_t unitary = eigensolver.eigenvectors();
        this->radial_ortho_matrix =
            unitary * eigs_invsqrt.matrix().asDiagonal() * unitary.adjoint();
      }

      std::shared_ptr<AtomicSmearingSpecificationBase> atomic_smearing{};
      AtomicSmearingType atomic_smearing_type{};
      math::Hyp1f1SphericalExpansion hyp1f1_calculator{true, 1e-13, 200};
      // data member used to store the contributions to the expansion
      Matrix_t radial_integral_neighbour{};
      Vector_t radial_integral_center{};

      Hypers_t hypers{};
      // some
      double interaction_cutoff{};
      size_t max_radial{};
      size_t max_angular{};

      // \sigma_n = (r_\text{cut}-\delta r_\text{cut})
      // \max(\sqrt{n},1)/n_\text{max}
      Vector_t radial_sigmas{};
      // b = 1 / (2*\sigma_n^2)
      Vector_t fac_b{};
      Vector_t radial_norm_factors{};
      Vector_t radial_n_factors{};
      Matrix_t radial_ortho_matrix{};
    };

  }  // namespace internal

  template <internal::RadialBasisType Type, class Hypers>
  decltype(auto) make_radial_integral(const Hypers & basis_hypers) {
    return std::static_pointer_cast<internal::RadialContributionBase>(
        std::make_shared<internal::RadialContribution<Type>>(basis_hypers));
  }

  template <internal::RadialBasisType Type>
  decltype(auto) downcast_radial_integral(
      const std::shared_ptr<internal::RadialContributionBase> &
          radial_integral) {
    return std::static_pointer_cast<internal::RadialContribution<Type>>(
        radial_integral);
  }

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
      using internal::RadialBasisType;

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

      // create the class that will compute the radial terms of the
      // expansion. the atomic smearing is an integral part of the
      // radial contribution
      if (radial_contribution_type.compare("GTO") == 0) {
        auto rc_shared = std::make_shared<
            internal::RadialContribution<RadialBasisType::GTO>>(hypers);
        this->atomic_smearing_type = rc_shared->atomic_smearing_type;
        this->radial_integral = rc_shared;
        this->radial_integral_type = RadialBasisType::GTO;

      } else {
        throw std::logic_error(
            "Requested Radial contribution type \'" + radial_contribution_type +
            "\' has not been implemented.  Must be one of" + ": \'GTO\'.");
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

    //! compute representation. choose the CutoffFunctionType from the hypers
    void compute();

    //! choose the RadialBasisType and AtomicSmearingType from the hypers
    template <internal::CutoffFunctionType FcType>
    void compute_by_radial_contribution();

    //! Compute the spherical exansion given several options
    template <internal::CutoffFunctionType FcType,
              internal::RadialBasisType RadialType,
              internal::AtomicSmearingType SmearingType>
    void compute_impl();

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

    internal::AtomicSmearingType atomic_smearing_type{};

    std::shared_ptr<internal::RadialContributionBase> radial_integral{};
    internal::RadialBasisType radial_integral_type{};

    std::shared_ptr<internal::CutoffFunctionBase> cutoff_function{};
    internal::CutoffFunctionType cutoff_function_type{};

    Hypers_t hypers{};
  };

  // compute classes template construction
  template <class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::compute() {
    // specialize based on the cutoff function
    using internal::CutoffFunctionType;

    switch (this->cutoff_function_type) {
    case CutoffFunctionType::Cosine: {
      this->compute_by_radial_contribution<CutoffFunctionType::Cosine>();
      break;
    }
    default:
      throw std::logic_error("The combination of parameter is not handdled.");
      break;
    }
  }

  template <class Mngr>
  template <internal::CutoffFunctionType FcType>
  void RepresentationManagerSphericalExpansion<
      Mngr>::compute_by_radial_contribution() {
    // specialize based on the type of radial contribution
    using internal::AtomicSmearingType;
    using internal::RadialBasisType;

    switch (internal::combineEnums(this->radial_integral_type,
                                   this->atomic_smearing_type)) {
    case internal::combineEnums(RadialBasisType::GTO,
                                AtomicSmearingType::Constant): {
      this->compute_impl<FcType, RadialBasisType::GTO,
                         AtomicSmearingType::Constant>();
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
  template <internal::CutoffFunctionType FcType,
            internal::RadialBasisType RadialType,
            internal::AtomicSmearingType SmearingType>
  void RepresentationManagerSphericalExpansion<Mngr>::compute_impl() {
    using math::PI;
    using math::pow;

    // downcast cutoff and radial contributions so they are functional
    auto cutoff_function{
        downcast_cutoff_function<FcType>(this->cutoff_function)};
    auto radial_integral{
        downcast_radial_integral<RadialType>(this->radial_integral)};

    auto n_row{this->max_radial};
    auto n_col{pow(this->max_angular + 1, 2)};
    this->expansions_coefficients.clear();
    this->expansions_coefficients.set_shape(n_row, n_col);
    this->expansions_coefficients.resize();

    for (auto center : this->structure_manager) {
      auto & coefficients_center = this->expansions_coefficients[center];
      Key_t center_type{center.get_atom_type()};

      //TODO(felix) think about an option to have "global" species,
      //"structure" species(or not), or automatic at the level of environment
      std::unordered_set<Key_t, internal::Hash<Key_t>> keys{};
      for (auto neigh : center) {
        keys.insert({neigh.get_atom_type()});
      }
      keys.insert({center_type});
      // initialize the expansion coefficients to 0
      coefficients_center.resize(keys, n_row, n_col);

      // Start the accumulator with the central atom
      coefficients_center[center_type].col(0) +=
          radial_integral->template compute_center_contribution<SmearingType>(
              center) /
          sqrt(4.0 * PI);

      for (auto neigh : center) {
        auto dist{this->structure_manager->get_distance(neigh)};
        auto direction{this->structure_manager->get_direction_vector(neigh)};
        Key_t neigh_type{neigh.get_atom_type()};

        // Note: the copy _should_ be optimized out (RVO)
        auto harmonics =
            math::compute_spherical_harmonics(direction, this->max_angular);
        auto neighbour_contribution =
            radial_integral
                ->template compute_neighbour_contribution<SmearingType>(dist,
                                                                        neigh);

        harmonics *= cutoff_function->f_c(dist);
        auto && coefficients_center_by_type{coefficients_center[neigh_type]};
        for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
          size_t lm_collective_idx{0};
          for (size_t angular_l{0}; angular_l < this->max_angular + 1;
               angular_l++) {
            for (size_t m_array_idx{0}; m_array_idx < 2 * angular_l + 1;
                 m_array_idx++) {
              coefficients_center_by_type(radial_n, lm_collective_idx) +=
                  //TODO(felix) check if you gain anything by accessing this
                  // outside the m loop
                  neighbour_contribution(radial_n, angular_l) *
                  harmonics(angular_l, m_array_idx);
                  //TODO(felix) refactor compute_spherical_harmonics to
                  //return a vector of (lmax+1)^2 elements
              ++lm_collective_idx;
            }
          }
        }
      }  // for (neigh : center)
    }    // for (center : structure_manager)
  }      // compute()

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_HH_
