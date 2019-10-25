/**
 * @file   calculator_spherical_expansion.hh
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

#ifndef SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_HH_
#define SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_HH_

#include "math/bessel.hh"
#include "math/gauss_legendre.hh"
#include "math/hyp1f1.hh"
#include "math/math_utils.hh"
#include "math/spherical_harmonics.hh"
#include "rascal_utility.hh"
#include "representations/calculator_base.hh"
#include "representations/cutoff_functions.hh"
#include "structure_managers/property_block_sparse.hh"
#include "structure_managers/structure_manager.hh"

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

  namespace internal {

    /**
     * List of possible Radial basis that can be used by the spherical
     * expansion.
     */
    enum class RadialBasisType { GTO, DVR, End_ };

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

      using Hypers_t = CalculatorBase::Hypers_t;
    };

    /**
     * Specification to hold the parameter for the atomic smearing function,
     * currently only Gaussians are supported.
     *
     * This is \f$\sigma\f$ in the definition
     * \f$f(r) = A \exp{\frac{-r^2}{2 \sigma^2}}\f$.
     * The width may depend both on the atomic species of the neighbour as well
     * as the distance.
     *
     * Note that this function is template-specialized by Gaussian sigma type
     * (constant, per-species, or radially dependent).
     *
     * @throw logic_error if the requested sigma type has not been implemented
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
    auto make_atomic_smearing(const Hypers & sigma_hypers) {
      return std::static_pointer_cast<AtomicSmearingSpecificationBase>(
          std::make_shared<AtomicSmearingSpecification<Type>>(sigma_hypers));
    }

    //! Utility to cast base to child class
    template <AtomicSmearingType Type>
    auto downcast_atomic_smearing(
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

      using Hypers_t = CalculatorBase::Hypers_t;
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
      // virtual Matrix_Ref compute_neighbour_derivative() = 0;
    };

    template <RadialBasisType RBT>
    struct RadialContribution {};

    /**
     * Implementation of the radial contribution for Gaussian Type Orbitals
     * radial basis functions and gaussian smearing of the atom density.
     *
     * @f[
     *      R^{GTO}_{n}(r) = \mathcal{N}_n\ r^{n} \exp[-br^2]
     * @f]
     *
     * @f{gather*}
     *      \newcommand{\dd}{\mathrm{d}\,}
     *      \mathcal{N}_n^2 = \frac{2}{\sigma_n^{2n + 3}\Gamma(n + 3/2)}\\
     *      \sigma_n = (r_\text{cut}-\delta r_\text{cut})
     *      \max(\sqrt{n},1)/n_\text{max}\\
     *      b=\frac{1}{2\sigma_n}\\
     *      \int_0^\infty R^{GTO}_{n}(r) R^{GTO}_{n^\prime}(r)
     *      \dd{r}= 2 \left(\frac{1}{2 \sigma_{n}^2}+
     *      \frac{1}{2 \sigma_{n^\prime}^2} \right)^{-\frac{1}{2}
     *      (3+n+n^\prime)} \Gamma(\frac{3+n+n^\prime}{2})
     * @f}
     *
     * See [the theory page](../SOAP.html#gto-like-radial-basis) for more
     * details.
     */
    template <>
    struct RadialContribution<RadialBasisType::GTO> : RadialContributionBase {
      //! Constructor
      explicit RadialContribution(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
        this->precompute();
      }
      // Destructor
      virtual ~RadialContribution() = default;
      // Copy constructor
      RadialContribution(const RadialContribution & other) = delete;
      // Move constructor
      RadialContribution(RadialContribution && other) = default;
      // Copy assignment operator
      RadialContribution & operator=(const RadialContribution & other) = delete;
      // Move assignment operator
      RadialContribution & operator=(RadialContribution && other) = default;

      using Parent = RadialContributionBase;
      using Hypers_t = typename Parent::Hypers_t;
      // using Matrix_t = typename Parent::Matrix_t;
      using Matrix_t = Eigen::MatrixXd;
      using Vector_t = typename Parent::Vector_t;
      using Matrix_Ref = typename Parent::Matrix_Ref;
      using Vector_Ref = typename Parent::Vector_Ref;

      /**
       * Set hyperparameters.
       * @param hypers is expected to be the same as the the input of
       *         the spherical expansion
       */
      void set_hyperparameters(const Hypers_t & hypers) {
        this->hypers = hypers;

        this->max_radial = hypers.at("max_radial");
        this->max_angular = hypers.at("max_angular");

        if (hypers.find("compute_gradients") != hypers.end()) {
          this->compute_gradients = hypers.at("compute_gradients").get<bool>();
        } else {  // Default false (don't compute gradients)
          this->compute_gradients = false;
        }

        // init size of the member data
        // both precomputed quantities and actual expansion coefficients
        this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
        this->ortho_norm_matrix.resize(this->max_radial, this->max_radial);
        this->fac_b.resize(this->max_radial, 1);
        this->a_b_l_n.resize(this->max_radial, this->max_angular + 1);
        this->distance_fac_a_l.resize(this->max_angular + 1);
        this->radial_norm_factors.resize(this->max_radial, 1);
        this->radial_n_factors.resize(this->max_radial);
        this->radial_sigmas.resize(this->max_radial, 1);
        this->radial_integral_neighbour.resize(this->max_radial,
                                               this->max_angular + 1);
        this->radial_integral_center.resize(this->max_radial);
        this->radial_neighbour_derivative.resize(this->max_radial,
                                                 this->max_angular + 1);

        // find the cutoff radius of the representation
        auto fc_hypers = hypers.at("cutoff_function").get<json>();
        this->interaction_cutoff =
            fc_hypers.at("cutoff").at("value").get<double>();

        // define the type of smearing to use
        auto smearing_hypers = hypers.at("gaussian_density").get<json>();
        auto smearing_type = smearing_hypers.at("type").get<std::string>();
        if (smearing_type == "Constant") {
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
        this->ortho_norm_matrix =
            this->radial_norm_factors.asDiagonal() * this->radial_ortho_matrix;

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
          double a_b_l_n{0};
          if (radial_n % 2 == 0) {
            a_b_l_n = sqrt(pow(fac_a + this->fac_b[radial_n],
                               -1 * static_cast<int>(3 + radial_n)));
          } else {
            a_b_l_n = pow(fac_a + this->fac_b[radial_n],
                          -1 * static_cast<int>(3 + radial_n) / 2);
          }

          this->radial_integral_center(radial_n) =
              this->radial_n_factors(radial_n) * a_b_l_n;
        }

        return Vector_Ref(this->radial_integral_center);
      }

      //! define the contribution from a neighbour atom to the expansion
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(double distance,
                                     ClusterRefKey<Order, Layer> & pair) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};
        // a = 1 / (2*\sigma^2)
        double fac_a{0.5 * pow(smearing->get_gaussian_sigma(pair), -2)};

        // computes (r_{ij}*a)^l incrementally
        this->distance_fac_a_l(0) = 1.;
        double distance_fac_a{distance * fac_a};
        for (size_t angular_l{1}; angular_l < this->max_angular + 1;
             angular_l++) {
          this->distance_fac_a_l(angular_l) =
              this->distance_fac_a_l(angular_l - 1) * distance_fac_a;
        }

        // computes (a+b_n)^{-0.5*(3+l+n)}
        Eigen::ArrayXd a_b_l{Eigen::rsqrt(fac_a + this->fac_b.array())};
        for (size_t radial_n{0}; radial_n < this->max_radial; radial_n++) {
          this->a_b_l_n(radial_n, 0) = pow(a_b_l(radial_n), 3 + radial_n);
        }
        // seems like vetorization does not improve things here because it is
        // memory is not contiguous ?
        for (size_t angular_l{1}; angular_l < this->max_angular + 1;
             ++angular_l) {
          this->a_b_l_n.col(angular_l) =
              (this->a_b_l_n.col(angular_l - 1).array() * a_b_l).matrix();
        }

        this->hyp1f1_calculator.calc(distance, fac_a, this->fac_b,
                                     this->compute_gradients);

        this->radial_integral_neighbour =
            (this->a_b_l_n.array() *
             this->hyp1f1_calculator.get_values().array())
                .matrix() *
            this->distance_fac_a_l.asDiagonal();

        return Matrix_Ref(this->radial_integral_neighbour);
      }

      /**
       * Compute the radial derivative of the neighbour contribution
       *
       * Note that you _must_ call compute_neighbour_contribution() first to
       * populate the relevant arrays!
       *
       * The derivative is taken with respect to the pair distance,
       * \f$r_{ij}\f$.  In order to get the radial component of the gradient,
       * remember to multiply by the direction vector
       * \f$
       *    \renewcommand{\vec}[1]{\mathbf{#1}}
       *    \hat{\vec{r}_{ij}}
       * \f$
       * (and not the vector itself), since
       * \f[
       *    \let\grad\nabla
       *    \grad_{\vec{r}_i} f(r_{ij}) =
       *                    \frac{\dd f}{\dd r_{ij}}
       *                    \frac{- \vec{r}_{ij}}{r_{ij}}
       *                  = \frac{\dd f}{\dd r_{ij}} -\hat{\vec{r}_{ij}}.
       * \f]
       *
       * so multiply by _negative_ \f$\hat{\vec{r}}_{ij}\f$ to get the radial
       * component of the gradient wrt motion of the central atom
       * (\f$\frac{d}{d\vec{r}_i}\f$).
       *
       * And finally, there is no compute_center_derivative() because that's
       * just zero -- the center contribution doesn't vary w.r.t. motion of
       * the central atom
       */
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_derivative(double distance,
                                   ClusterRefKey<Order, Layer> & /*pair*/) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        // start proportional_factors as the list of l from 0 to l_max
        Vector_t proportional_factors =
            Vector_t::LinSpaced(this->max_angular + 1, 0, this->max_angular);
        proportional_factors /= distance;

        this->radial_neighbour_derivative =
            (this->a_b_l_n.array() *
             this->hyp1f1_calculator.get_derivatives().array())
                .matrix() *
            this->distance_fac_a_l.asDiagonal();

        this->radial_neighbour_derivative +=
            this->radial_integral_neighbour * proportional_factors.asDiagonal();

        return Matrix_Ref(this->radial_neighbour_derivative);
      }

      template <typename Coeffs>
      void finalize_coefficients(Coeffs & coefficients) const {
        coefficients.lhs_dot(this->ortho_norm_matrix);
      }

      template <int n_spatial_dimensions, typename Coeffs, typename Center>
      void finalize_coefficients_der(Coeffs & coefficients_gradient,
                                     Center & center) const {
        for (auto neigh : center.with_self_pair()) {
          auto & coefficients_neigh_gradient = coefficients_gradient[neigh];
          coefficients_neigh_gradient
              .template lhs_dot_der<n_spatial_dimensions>(
                  this->ortho_norm_matrix);
        }  // for (neigh : center)
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
        // UPDATE nevermind, the overlap matrix becomes singular well before
        // this point.
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
        Eigen::ArrayXd eigs_invsqrt = eigenvalues.array().rsqrt();
        Matrix_t unitary = eigensolver.eigenvectors();
        this->radial_ortho_matrix =
            unitary * eigs_invsqrt.matrix().asDiagonal() * unitary.adjoint();
      }

      Matrix_t get_radial_orthonormalization_matrix() const {
        return this->radial_norm_factors.asDiagonal() *
               this->radial_ortho_matrix;
      }

      Matrix_Ref get_radial_integral_neighbour() const {
        return Matrix_Ref(this->radial_integral_neighbour);
      }

      Matrix_Ref get_radial_neighbour_derivative() const {
        return Matrix_Ref(this->radial_neighbour_derivative);
      }

      std::shared_ptr<AtomicSmearingSpecificationBase> atomic_smearing{};
      AtomicSmearingType atomic_smearing_type{};
      math::Hyp1f1SphericalExpansion hyp1f1_calculator{true, 1e-13, 200};
      // data member used to store the contributions to the expansion
      Matrix_t radial_integral_neighbour{};
      Vector_t radial_integral_center{};
      // And derivatives
      Matrix_t radial_neighbour_derivative{};
      // and of course, d/dr of the center contribution is zero

      Hypers_t hypers{};
      // some usefull parameters
      double interaction_cutoff{};
      size_t max_radial{};
      size_t max_angular{};
      bool compute_gradients{};

      // \sigma_n = (r_\text{cut}-\delta r_\text{cut})
      // \max(\sqrt{n},1)/n_\text{max}
      Vector_t radial_sigmas{};
      // b = 1 / (2*\sigma_n^2)
      Vector_t fac_b{};
      Matrix_t a_b_l_n{};
      Vector_t distance_fac_a_l{};
      Vector_t radial_norm_factors{};
      Vector_t radial_n_factors{};
      Matrix_t radial_ortho_matrix{};
      Matrix_t ortho_norm_matrix{};
    };

    /**
     * Implementation of the radial contribution for DVR basis
     *
     * See the
     * [theory page](../SOAP.html#numerical-integration-of-the-radial-integral)
     * for more details.
     */
    template <>
    struct RadialContribution<RadialBasisType::DVR> : RadialContributionBase {
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
       *
       * @param hypers is expected to be the same as the the input of
       *         the spherical expansion
       */
      void set_hyperparameters(const Hypers_t & hypers) {
        this->hypers = hypers;

        this->max_radial = hypers.at("max_radial");
        this->max_angular = hypers.at("max_angular");

        // init size of the member data
        // both precomputed quantities and actual expansion coefficients
        this->legendre_radial_factor.resize(this->max_radial);
        this->legendre_points.resize(this->max_radial);

        this->radial_integral_neighbour.resize(this->max_radial,
                                               this->max_angular + 1);
        this->radial_neighbour_derivative.resize(this->max_radial,
                                                 this->max_angular + 1);
        this->radial_integral_center.resize(this->max_radial);

        // find the cutoff radius of the representation
        auto fc_hypers = hypers.at("cutoff_function").get<json>();
        this->interaction_cutoff =
            fc_hypers.at("cutoff").at("value").get<double>();
        this->smooth_width =
            fc_hypers.at("smooth_width").at("value").get<double>();

        // define the type of smearing to use
        auto smearing_hypers = hypers.at("gaussian_density").get<json>();
        auto smearing_type = smearing_hypers.at("type").get<std::string>();
        if (smearing_type == "Constant") {
          this->atomic_smearing_type = AtomicSmearingType::Constant;
          this->atomic_smearing =
              make_atomic_smearing<AtomicSmearingType::Constant>(
                  smearing_hypers);
          this->smearing =
              smearing_hypers.at("gaussian_sigma").at("value").get<double>();
        } else {
          throw std::logic_error(
              "Requested Gaussian sigma type \'" + smearing_type +
              "\' has not been implemented.  Must be one of" +
              ": \'Constant\'.");
        }
      }

      void precompute() {
        auto point_weight{math::compute_gauss_legendre_points_weights(
            0., this->interaction_cutoff + 3 * this->smearing,
            this->max_radial)};

        // sqrt(w) * x
        // (if you think it should be x^2 and not x, think again -- the
        // transformation from integrating the overlap in 3-D spherical
        // coordinates to the 1-D radial coordinate absorbs a factor of r.
        // For more details, see the SOAP theory documentation)
        // TODO(max) link to SOAP theory documentation
        this->legendre_radial_factor =
            point_weight.col(1).array().sqrt() * point_weight.col(0).array();

        this->legendre_points = point_weight.col(0);

        this->bessel.precompute(this->max_angular, this->legendre_points);
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

        this->radial_integral_center =
            this->legendre_radial_factor.array() *
            Eigen::exp((-fac_a * this->legendre_points.array().square()));

        return Matrix_Ref(this->radial_integral_center);
      }

      //! define the contribution from a neighbour atom to the expansion
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(double distance,
                                     ClusterRefKey<Order, Layer> & pair) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};
        // a = 1 / (2*\sigma^2)
        double fac_a{0.5 * pow(smearing->get_gaussian_sigma(pair), -2)};

        this->bessel.calc(distance, fac_a);

        this->radial_integral_neighbour =
            this->legendre_radial_factor.asDiagonal() *
            this->bessel.get_values().matrix();

        return Matrix_Ref(this->radial_integral_neighbour);
      }

      /**
       * Compute the radial derivative of the neighbour contribution
       *
       * @todo still needs to be implemented for the DVR radial basis
       */
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_derivative(double /*distance*/,
                                   ClusterRefKey<Order, Layer> & /*pair*/) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        // TODO(max,felix) implement (!) -- these are dummy values
        return Matrix_Ref(this->radial_neighbour_derivative);
      }

      template <typename Coeffs>
      void finalize_coefficients(Coeffs & /*coefficients*/) const {}

      template <int n_spatial_dimensions, typename Coeffs, typename Center>
      void finalize_coefficients_der(Coeffs & /*coefficients_gradient*/,
                                     Center & /*center*/) const {}

      math::ModifiedSphericalBessel bessel{};

      std::shared_ptr<AtomicSmearingSpecificationBase> atomic_smearing{};
      AtomicSmearingType atomic_smearing_type{};

      // data member used to store the contributions to the expansion
      Matrix_t radial_integral_neighbour{};
      Matrix_t radial_neighbour_derivative{};
      Vector_t radial_integral_center{};

      Hypers_t hypers{};
      // some useful parameters
      double interaction_cutoff{};
      double smooth_width{};
      double smearing{};
      size_t max_radial{};
      size_t max_angular{};

      Vector_t legendre_radial_factor{};
      Vector_t legendre_points{};
      Vector_t legendre_points2{};
    };

  }  // namespace internal

  template <internal::RadialBasisType Type, class Hypers>
  auto make_radial_integral(const Hypers & basis_hypers) {
    return std::static_pointer_cast<internal::RadialContributionBase>(
        std::make_shared<internal::RadialContribution<Type>>(basis_hypers));
  }

  template <internal::RadialBasisType Type>
  auto downcast_radial_integral(
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
   * of spherical harmonics (à la SOAP) and a radial basis of
   * either Gaussians (again, as in SOAP/SphericalInvariants) or one of the more
   * recent bases currently under development.
   */
  class CalculatorSphericalExpansion : public CalculatorBase {
   public:
    using Parent = CalculatorBase;
    using Hypers_t = typename Parent::Hypers_t;
    using ReferenceHypers_t = Parent::ReferenceHypers_t;
    using Key_t = typename Parent::Key_t;

    template <class StructureManager>
    using Property_t =
        BlockSparseProperty<double, 1, 0, StructureManager, Key_t>;
    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, 0, StructureManager, Key_t>;

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
    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::AtomicSmearingType;
      using internal::CutoffFunctionType;
      using internal::RadialBasisType;
      this->hypers = hypers;

      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      if (hypers.find("n_species") != hypers.end()) {
        this->n_species = hypers.at("n_species");
      } else {
        this->n_species = 1;  // default: no species distinction
      }
      if (hypers.find("compute_gradients") != hypers.end()) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }
      this->spherical_harmonics.precompute(this->max_angular,
                                           this->compute_gradients);

      auto radial_contribution_hypers =
          hypers.at("radial_contribution").get<json>();
      auto radial_contribution_type =
          radial_contribution_hypers.at("type").get<std::string>();

      // create the class that will compute the radial terms of the
      // expansion. the atomic smearing is an integral part of the
      // radial contribution
      if (radial_contribution_type == "GTO") {
        auto rc_shared = std::make_shared<
            internal::RadialContribution<RadialBasisType::GTO>>(hypers);
        this->atomic_smearing_type = rc_shared->atomic_smearing_type;
        this->radial_integral = rc_shared;
        this->radial_integral_type = RadialBasisType::GTO;

      } else if (radial_contribution_type == "DVR") {
        auto rc_shared = std::make_shared<
            internal::RadialContribution<RadialBasisType::DVR>>(hypers);
        this->atomic_smearing_type = rc_shared->atomic_smearing_type;
        this->radial_integral = rc_shared;
        this->radial_integral_type = RadialBasisType::DVR;
      } else {
        throw std::logic_error("Requested Radial contribution type \'" +
                               radial_contribution_type +
                               "\' has not been implemented.  Must be one of" +
                               ": \'GTO\' or \'DVR\'. ");
      }

      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      auto fc_type = fc_hypers.at("type").get<std::string>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
      this->cutoff_smooth_width = fc_hypers.at("smooth_width").at("value");
      if (fc_type == "ShiftedCosine") {
        this->cutoff_function_type = CutoffFunctionType::ShiftedCosine;
        this->cutoff_function =
            make_cutoff_function<CutoffFunctionType::ShiftedCosine>(fc_hypers);
      } else if (fc_type == "RadialScaling") {
        this->cutoff_function_type = CutoffFunctionType::RadialScaling;
        this->cutoff_function =
            make_cutoff_function<CutoffFunctionType::RadialScaling>(fc_hypers);
      } else {
        throw std::logic_error("Requested cutoff function type \'" + fc_type +
                               "\' has not been implemented.  Must be one of" +
                               ": \'ShiftedCosine\' or 'RadialScaling'.");
      }

      this->set_name(hypers);
    }

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
    explicit CalculatorSphericalExpansion(const Hypers_t & hyper)
        : CalculatorBase{} {
      this->set_default_prefix("spherical_expansion_");
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    CalculatorSphericalExpansion(const CalculatorSphericalExpansion & other) =
        delete;

    //! Move constructor
    CalculatorSphericalExpansion(CalculatorSphericalExpansion && other) =
        default;

    //! Destructor
    virtual ~CalculatorSphericalExpansion() = default;

    //! Copy assignment operator
    CalculatorSphericalExpansion &
    operator=(const CalculatorSphericalExpansion & other) = delete;

    //! Move assignment operator
    CalculatorSphericalExpansion &
    operator=(CalculatorSphericalExpansion && other) = default;

    /**
     * Compute representation for a given structure manager.
     *
     * @tparam StructureManager a (single or collection)
     * of structure manager(s) (in an iterator) held in shared_ptr
     */
    template <class StructureManager>
    void compute(StructureManager & managers);

    //! choose the RadialBasisType and AtomicSmearingType from the hypers
    template <internal::CutoffFunctionType FcType, class StructureManager>
    void compute_by_radial_contribution(StructureManager & managers);

    /**
     * loop over a collection of manangers if it is an iterator.
     * Or just call compute_impl() if it's a single manager (see below)
     */
    template <
        internal::CutoffFunctionType FcType,
        internal::RadialBasisType RadialType,
        internal::AtomicSmearingType SmearingType, class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<FcType, RadialType, SmearingType>(manager);
      }
    }

    //! single manager case
    template <internal::CutoffFunctionType FcType,
              internal::RadialBasisType RadialType,
              internal::AtomicSmearingType SmearingType, class StructureManager,
              std::enable_if_t<
                  not(internal::is_proper_iterator<StructureManager>::value),
                  int> = 0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<FcType, RadialType, SmearingType>(manager);
    }

    //! Compute the spherical exansion given several options
    template <internal::CutoffFunctionType FcType,
              internal::RadialBasisType RadialType,
              internal::AtomicSmearingType SmearingType, class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

   protected:
    double interaction_cutoff{};
    double cutoff_smooth_width{};
    size_t max_radial{};
    size_t max_angular{};
    size_t n_species{};
    bool compute_gradients{};

    internal::AtomicSmearingType atomic_smearing_type{};

    std::shared_ptr<internal::RadialContributionBase> radial_integral{};
    internal::RadialBasisType radial_integral_type{};

    std::shared_ptr<internal::CutoffFunctionBase> cutoff_function{};
    internal::CutoffFunctionType cutoff_function_type{};

    Hypers_t hypers{};

    math::SphericalHarmonics spherical_harmonics{};
  };

  // compute classes template construction
  template <class StructureManager>
  void CalculatorSphericalExpansion::compute(StructureManager & managers) {
    // specialize based on the cutoff function
    using internal::CutoffFunctionType;

    switch (this->cutoff_function_type) {
    case CutoffFunctionType::ShiftedCosine: {
      this->compute_by_radial_contribution<CutoffFunctionType::ShiftedCosine>(
          managers);
      break;
    }
    case CutoffFunctionType::RadialScaling: {
      this->compute_by_radial_contribution<CutoffFunctionType::RadialScaling>(
          managers);
      break;
    }
    default:
      // The control flow really should never reach here.  But just in case,
      // provide the necessary information to debug this problem.
      std::basic_ostringstream<char> err_message;
      err_message << "Invalid cutoff function type encountered ";
      err_message << "(This is a bug.  Debug info for developers: ";
      err_message << "cutoff_function_type == ";
      err_message << static_cast<int>(this->cutoff_function_type);
      err_message << ")" << std::endl;
      throw std::logic_error(err_message.str());
      break;
    }
  }

  template <internal::CutoffFunctionType FcType, class StructureManager>
  void CalculatorSphericalExpansion::compute_by_radial_contribution(
      StructureManager & managers) {
    // specialize based on the type of radial contribution
    using internal::AtomicSmearingType;
    using internal::RadialBasisType;

    switch (internal::combineEnums(this->radial_integral_type,
                                   this->atomic_smearing_type)) {
    case internal::combineEnums(RadialBasisType::GTO,
                                AtomicSmearingType::Constant): {
      this->compute_loop<FcType, RadialBasisType::GTO,
                         AtomicSmearingType::Constant>(managers);
      break;
    }
    case internal::combineEnums(RadialBasisType::DVR,
                                AtomicSmearingType::Constant): {
      this->compute_loop<FcType, RadialBasisType::DVR,
                         AtomicSmearingType::Constant>(managers);
      break;
    }
    default:
      // The control flow really should never reach here.  In this case, any
      // "invalid combination of parameters" should have already been handled at
      // the parameter processing stage where the user can be notified in a
      // helpful way.  But in case we do get here, provide the necessary
      // information to debug this problem.
      std::basic_ostringstream<char> err_message;
      err_message << "Invalid combination of atomic smearing and radial basis ";
      err_message << "type encountered (This is a bug.  Debug info for ";
      err_message << "developers: "
                  << "radial_integral_type == ";
      err_message << static_cast<int>(this->radial_integral_type);
      err_message << ", atomic_smearing_type == ";
      err_message << static_cast<int>(this->atomic_smearing_type);
      err_message << ")" << std::endl;
      throw std::logic_error(err_message.str());
      break;
    }
  }

  /**
   * Compute the spherical expansion
   */
  template <internal::CutoffFunctionType FcType,
            internal::RadialBasisType RadialType,
            internal::AtomicSmearingType SmearingType, class StructureManager>
  void CalculatorSphericalExpansion::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using Prop_t = Property_t<StructureManager>;
    using PropGrad_t = PropertyGradient_t<StructureManager>;
    constexpr static int n_spatial_dimensions = StructureManager::dim();

    using math::PI;
    using math::pow;

    auto && expansions_coefficients{
        manager->template get_property_ref<Prop_t>(this->get_name())};

    auto && expansions_coefficients_gradient{
        manager->template get_property_ref<PropGrad_t>(
            this->get_gradient_name())};

    // if the representation has already been computed for the current
    // structure then do nothing
    if (expansions_coefficients.is_updated()) {
      return;
    }

    // downcast cutoff and radial contributions so they are functional
    auto cutoff_function{
        downcast_cutoff_function<FcType>(this->cutoff_function)};
    auto radial_integral{
        downcast_radial_integral<RadialType>(this->radial_integral)};

    auto n_row{this->max_radial};
    auto n_col{(this->max_angular + 1) * (this->max_angular + 1)};
    expansions_coefficients.clear();
    expansions_coefficients.set_shape(n_row, n_col);
    expansions_coefficients.resize();

    if (this->compute_gradients) {
      expansions_coefficients_gradient.clear();
      // Row-major ordering, so the Cartesian (spatial) index varies slowest
      expansions_coefficients_gradient.set_shape(n_spatial_dimensions * n_row,
                                                 n_col);
      expansions_coefficients_gradient.resize();
    }

    /* @TODO(felix,max) use the parity of the spherical harmonics to use half
     * neighbourlist, i.e. C^{ij}_{nlm} = (-1)^l C^{ji}_{nlm}.
     */
    for (auto center : manager) {
      auto & coefficients_center = expansions_coefficients[center];
      auto & coefficients_center_gradient =
          expansions_coefficients_gradient[center.get_atom_ii()];
      Key_t center_type{center.get_atom_type()};

      // TODO(felix) think about an option to have "global" species,
      // "structure" species(or not), or automatic at the level of environment
      std::unordered_set<Key_t, internal::Hash<Key_t>> keys{};
      for (auto neigh : center) {
        keys.insert({neigh.get_atom_type()});
      }
      keys.insert({center_type});
      // initialize the expansion coefficients to 0
      coefficients_center.resize(keys, n_row, n_col, 0.);
      if (this->compute_gradients) {
        coefficients_center_gradient.resize(keys, n_spatial_dimensions * n_row,
                                            n_col, 0.);
      }

      // Start the accumulator with the central atom
      coefficients_center[center_type].col(0) +=
          radial_integral->template compute_center_contribution<SmearingType>(
              center) /
          sqrt(4.0 * PI);

      for (auto neigh : center) {
        auto dist{manager->get_distance(neigh)};
        auto direction{manager->get_direction_vector(neigh)};
        Key_t neigh_type{neigh.get_atom_type()};
        auto & coefficients_neigh_gradient =
            expansions_coefficients_gradient[neigh];

        this->spherical_harmonics.calc(direction, this->compute_gradients);
        auto && harmonics{spherical_harmonics.get_harmonics()};
        auto && harmonics_gradients{
            spherical_harmonics.get_harmonics_derivatives()};

        auto && neighbour_contribution =
            radial_integral
                ->template compute_neighbour_contribution<SmearingType>(dist,
                                                                        neigh);
        double f_c{cutoff_function->f_c(dist)};
        auto && coefficients_center_by_type{coefficients_center[neigh_type]};

        // compute the coefficients
        size_t l_block_idx{0};
        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             ++angular_l) {
          size_t l_block_size{2 * angular_l + 1};
          coefficients_center_by_type.block(0, l_block_idx, max_radial,
                                            l_block_size) +=
              (neighbour_contribution.col(angular_l) *
               (harmonics.segment(l_block_idx, l_block_size) * f_c));
          l_block_idx += l_block_size;
        }

        // compute the gradients of the coefficients with respect to
        // atoms positions
        if (this->compute_gradients) {
          std::vector<Key_t> neigh_types{neigh_type};
          coefficients_neigh_gradient.resize(
              neigh_types, n_spatial_dimensions * n_row, n_col, 0.);

          auto && neighbour_derivative =
              radial_integral
                  ->template compute_neighbour_derivative<SmearingType>(dist,
                                                                        neigh);
          double df_c{cutoff_function->df_c(dist)};
          // The gradients only contribute to the type of the neighbour
          // (the atom that's moving)
          // grad_i c^{ij}
          auto && gradient_center_by_type{
              coefficients_center_gradient[neigh_type]};
          // grad_j c^{ij}
          auto && gradient_neigh_by_type{
              coefficients_neigh_gradient[neigh_type]};

          // Radial component: d/dr_{ij} (c_{ij} f_c{r_{ij}}) \hat{r_{ij}}
          // clang-format off
          Matrix_t pair_gradient_contribution_p1 =
                ((neighbour_derivative * f_c)
                 + (neighbour_contribution * df_c));
          Matrix_t pair_gradient_contribution{this->max_radial,
                                              this->max_angular + 1};
          for (int cartesian_idx{0}; cartesian_idx < n_spatial_dimensions;
                 ++cartesian_idx) {
            size_t l_block_idx{0};
            for (size_t angular_l{0}; angular_l < this->max_angular + 1;
                ++angular_l) {
              size_t l_block_size{2 * angular_l + 1};
              pair_gradient_contribution.resize(this->max_radial, l_block_size);
              pair_gradient_contribution =
                pair_gradient_contribution_p1.col(angular_l)
                * harmonics.segment(l_block_idx, l_block_size)
                * direction(cartesian_idx);
              pair_gradient_contribution +=
                  neighbour_contribution.col(angular_l)
                  * harmonics_gradients.block(cartesian_idx, l_block_idx,
                                              1, l_block_size)
                  * f_c / dist;

              // Each Cartesian gradient component occupies a contiguous block
              // (row-major storage)
              gradient_center_by_type.block(
                  cartesian_idx * max_radial, l_block_idx,
                  max_radial, l_block_size) -= pair_gradient_contribution;
              gradient_neigh_by_type.block(
                  cartesian_idx * max_radial, l_block_idx,
                  max_radial, l_block_size) += pair_gradient_contribution;
              l_block_idx += l_block_size;
              // clang-format on
            }  // for (angular_l)
          }    // for cartesian_idx
        }      // if (this->compute_gradients)
      }        // for (neigh : center)

      // Normalize and orthogonalize the radial coefficients
      radial_integral->finalize_coefficients(coefficients_center);
      if (this->compute_gradients) {
        radial_integral
            ->template finalize_coefficients_der<n_spatial_dimensions>(
                expansions_coefficients_gradient, center);
      }
    }  // for (center : manager)
  }    // compute()

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_EXPANSION_HH_
