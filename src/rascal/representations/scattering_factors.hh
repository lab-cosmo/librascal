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

#ifndef SRC_RASCAL_REPRESENTATIONS_GEOMETRIC_KSPACE_FACTORS_HH_
#define SRC_RASCAL_REPRESENTATIONS_GEOMETRIC_KSPACE_FACTORS_HH_

#include "rascal/math/hyp1f1.hh"
#include "rascal/math/interpolator.hh"
#include "rascal/math/spherical_harmonics.hh"
#include "rascal/math/kvec_generator.hh"
#include "rascal/math/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/cutoff_functions.hh"
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

  namespace internal {

    /**
     * Base class to define the radial contribution to the spherical expansion
     */
    struct RadialContributionKspaceBase {
      //! Constructor
      RadialContributionKspaceBase() = default;
      //! Destructor
      virtual ~RadialContributionKspaceBase() = default;
      //! Copy constructor
      RadialContributionKspaceBase(const RadialContributionKspaceBase & other) = delete;
      //! Move constructor
      RadialContributionKspaceBase(RadialContributionKspaceBase && other) = default;
      //! Copy assignment operator
      RadialContributionKspaceBase &
      operator=(const RadialContributionKspaceBase & other) = delete;
      //! Move assignment operator
      RadialContributionKspaceBase &
      operator=(RadialContributionKspaceBase && other) = default;

      using Hypers_t = CalculatorBase::Hypers_t;
      using Matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;
      using Vector_t = Eigen::VectorXd;
      using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;
      using Vector_Ref = typename Eigen::Ref<const Vector_t>;

      //! Pure Virtual Function to set hyperparameters of the cutoff function
      virtual void set_hyperparameters(const Hypers_t &) = 0;

      virtual void precompute() = 0;
    };

    template <RadialBasisType RBT>
    struct RadialContributionKspace {};

    /**
     * Implementation of the radial contribution for Gaussian Type Orbitals
     * radial basis functions and gaussian smearing of the atom density.
     */
    template <>
    struct RadialContributionKspace<RadialBasisType::GTO> : RadialContributionKspaceBase {
      // Constructor
      explicit RadialContributionKspace(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
        this->precompute();
      }
      // Destructor
      virtual ~RadialContributionKspace() = default;
      // Copy constructor
      RadialContributionKspace(const RadialContributionKspace & other) = delete;
      // Move constructor
      RadialContributionKspace(RadialContributionKspace && other) = default;
      // Copy assignment operator
      RadialContributionKspace & operator=(const RadialContributionKspace & other) = delete;
      // Move assignment operator
      RadialContributionKspace & operator=(RadialContributionKspace && other) = default;

      using Parent = RadialContributionKspaceBase;
      using Hypers_t = typename Parent::Hypers_t;
      using Matrix_t = Eigen::MatrixXd;
      using Vector_t = typename Parent::Vector_t;
      using Matrix_Ref = typename Parent::Matrix_Ref;
      using Vector_Ref = typename Parent::Vector_Ref;

      /**
       * Set hyperparameters by overriding them from input ones.
       */
      void set_hyperparameters(const Hypers_t & hypers) override {
	this->hypers = hypers;
        this->max_radial = hypers.at("max_radial");
        this->max_angular = hypers.at("max_angular");
        // Gradient?
        if (hypers.find("compute_gradients") != hypers.end()) {
          this->compute_gradients = hypers.at("compute_gradients").get<bool>();
        } else {  
          this->compute_gradients = false;
        }

        // Initialize size of the member arrays
        this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
        this->ortho_norm_matrix.resize(this->max_radial, this->max_radial);
        this->radial_norm_factors.resize(this->max_radial, 1);
        this->radial_sigmas.resize(this->max_radial, 1);
        this->radial_nl_factors.resize(this->max_radial, this->max_angular + 1);
        this->radial_integral_neighbour.resize(this->max_radial,
                                               this->max_angular + 1);

        // find the cutoff radius of the representation
        auto fc_hypers = hypers.at("cutoff_function").get<json>();
        this->interaction_cutoff = fc_hypers.at("cutoff").at("value").get<double>();

        // define the type of smearing to use
        auto smearing_hypers = hypers.at("gaussian_density").get<json>();
        auto smearing_type = smearing_hypers.at("type").get<std::string>();
        if (smearing_type == "Constant") {
          this->atomic_smearing_type = AtomicSmearingType::Constant;
          this->atomic_smearing = make_atomic_smearing<AtomicSmearingType::Constant>(smearing_hypers);
        } else {
          throw std::logic_error(
              "Requested Gaussian sigma type \'" + smearing_type +
              "\' has not been implemented.  Must be one of: \'Constant\'.");
        }
      }

      void precompute() override {
        this->precompute_radial_prefactors();
        this->precompute_radial_overlap();
        this->ortho_norm_matrix = this->radial_norm_factors.asDiagonal() * this->radial_ortho_matrix;
        this->hyp1f1_calculator.precompute(this->max_radial, this->max_angular);
      }

      /**
       * Define the contribution from a neighbour atom to the expansion
       * without requiring a cluster object so it can be used with the interpolator.
       */
      template <AtomicSmearingType AST>
      Matrix_t compute_contribution(const double distance, const double sigma) {
        using math::PI;
        using math::pow;
        using std::sqrt;

        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};

        this->hyp1f1_calculator.calc(distance, fac_a, this->fac_b);

        Matrix_t radial_integral_neighbour =
            (this->hyp1f1_calculator.get_values().array())
	    .matrix() ;
        radial_integral_neighbour.transpose() *= this->radial_norm_factors.asDiagonal();
        return radial_integral_neighbour;
      }

      //! define the contribution from a neighbour atom to the expansion
      template <AtomicSmearingType AST, size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(const double distance,
                                     const ClusterRefKey<Order, Layer> & pair) {
        auto smearing{downcast_atomic_smearing<AST>(this->atomic_smearing)};
        double smearing_value{smearing->get_gaussian_sigma(pair)};
        return this->compute_neighbour_contribution(distance, fac_a);
      }

      // Define the contribution from a neighbour atom to the expansion with a precomputed a-factor
      Matrix_Ref compute_neighbour_contribution(const double distance,
                                                const double fac_a) {
        using math::pow;
        using std::sqrt;

        this->hyp1f1_calculator.calc(distance, fac_a, this->fac_b,
                                     this->compute_gradients);

        this->radial_integral_neighbour =
            (this->hyp1f1_calculator.get_values().array())
                .matrix() ;

        return Matrix_Ref(this->radial_integral_neighbour);
      }

      template <typename Coeffs>
      void finalize_coefficients(Coeffs & coefficients) const {
        coefficients.lhs_dot(this->ortho_norm_matrix);
      }

      template <int NDims, typename Coeffs, typename Center>
      void finalize_coefficients_der(Coeffs & coefficients_gradient,
                                     Center & center) const {
        for (auto neigh : center.pairs_with_self_pair()) {
          auto & coefficients_neigh_gradient = coefficients_gradient[neigh];
          coefficients_neigh_gradient.template lhs_dot_der<NDims>(
              this->ortho_norm_matrix);
        }  // for (neigh : center)
      }

      /** Compute common prefactors for the radial Gaussian basis functions */
      void precompute_radial_prefactors() {
        using math::pow;
        // Precompute common prefactors
        for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          this->radial_sigmas[radial_n] =
              std::max(std::sqrt(static_cast<double>(radial_n)), 1.0) *
              this->interaction_cutoff / static_cast<double>(this->max_radial);
          this->radial_norm_factors(radial_n) =
              std::sqrt(2.0 / (std::tgamma(1.5 + radial_n) *
                        pow(this->radial_sigmas[radial_n], 3 + 2 * radial_n)));
          for (size_t angular_l{0}; angular_l < this->max_angular+1; ++angular_l) {
             this->radial_nl_factors(radial_n,angular_l) =
		 pow(2.0,0.5*(radial_n - angular_l - 1)) 
                 * pow(this->radial_sigmas[radial_n], 3 + radial_n + angular_l)
                 * std::tgamma(0.5 * (3.0 + radial_n + angular_l)) / std::tgamma(1.5+angular_l);
          }
        }
      }

      /**
       * Compute the radial overlap matrix for later orthogonalization.
       * @throw runtime_error if the overlap matrix cannot be diagonalized.
       */
      void precompute_radial_overlap() {
        using math::pow;
        using std::sqrt;
        using std::tgamma;

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

      std::shared_ptr<AtomicSmearingSpecificationBase> atomic_smearing{};
      AtomicSmearingType atomic_smearing_type{};
      math::Hyp1f1SphericalExpansion hyp1f1_calculator{true, 1e-13, 200};
      // data member used to store the contributions to the expansion
      Matrix_t radial_integral_neighbour{};

      Hypers_t hypers{};
      // some usefull parameters
      double interaction_cutoff{};
      size_t max_radial{};
      size_t max_angular{};
      bool compute_gradients{};

      // \sigma_n = (r_\text{cut}-\delta r_\text{cut})
      // \max(\sqrt{n},1)/n_\text{max}
      Vector_t radial_sigmas{};
      Vector_t radial_norm_factors{};
      Vector_t radial_nl_factors{};
      Matrix_t radial_ortho_matrix{};
      Matrix_t ortho_norm_matrix{};
    };

    /* A RadialContributionKspaceHandler handles the different cases of
     * AtomicSmearingType and OptimizationType. Depending on these template
     * parameters different member variables have to be used and different
     * parameters can be precomputed.
     */
    template <RadialBasisType RBT, AtomicSmearingType AST, OptimizationType IT>
    struct RadialContributionKspaceHandler {};

    /* For the a constant smearing type the "a" factor can be precomputed
     */
    template <RadialBasisType RBT>
    struct RadialContributionKspaceHandler<RBT, AtomicSmearingType::Constant,
                                     OptimizationType::None>
        : public RadialContributionKspace<RBT> {
     public:
      using Parent = RadialContributionKspace<RBT>;
      using Hypers_t = typename Parent::Hypers_t;
      using Matrix_Ref = typename Parent::Matrix_Ref;
      using Vector_Ref = typename Parent::Vector_Ref;

      explicit RadialContributionKspaceHandler(const Hypers_t & hypers)
          : Parent(hypers) {
        this->precompute_fac_a();
      }

      template <size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(const double distance,
                                     const ClusterRefKey<Order, Layer> &) {
        return Parent::compute_neighbour_contribution(distance, this->fac_a);
      }

      Matrix_Ref compute_neighbour_contribution(const double) {
        return Parent::compute_neighbour_contribution(this->fac_a);
      }

      template <size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_derivative(const double distance,
                                   const ClusterRefKey<Order, Layer> & pair) {
        return Parent::compute_neighbour_derivative(distance, pair);
      }

     protected:
      void precompute_fac_a() {
        auto smearing{downcast_atomic_smearing<AtomicSmearingType::Constant>(
            this->atomic_smearing)};
        this->fac_a = 0.5 * pow(smearing->get_gaussian_sigma(), -2);
      }

      // 1/(2σ^2)
      double fac_a{};
    };

    /* For the a constant smearing type the "a" factor can be precomputed and
     * when using the interpolator has to be initialized and used.
     */
    template <RadialBasisType RBT>
    struct RadialContributionKspaceHandler<RBT, AtomicSmearingType::Constant,
                                     OptimizationType::Interpolator>
        : public RadialContributionKspace<RBT> {
     public:
      using Parent = RadialContributionKspace<RBT>;
      using Hypers_t = typename Parent::Hypers_t;
      using Matrix_t = typename Parent::Matrix_t;
      using Matrix_Ref = typename Parent::Matrix_Ref;
      using Vector_Ref = typename Parent::Vector_Ref;
      using Interpolator_t = math::InterpolatorMatrixUniformCubicSpline<
          math::RefinementMethod_t::Exponential>;

      explicit RadialContributionKspaceHandler(const Hypers_t & hypers)
          : Parent(hypers) {
        this->precompute();
        this->init_interpolator(hypers);
      }

      // If we find a case where smarter parameters for x1 and x2 can be given
      explicit RadialContributionKspaceHandler(const Hypers_t & hypers,
                                         const double x1, const double x2,
                                         const double accuracy)
          : Parent(hypers) {
        this->precompute();
        this->init_interpolator(x1, x2, accuracy);
      }

      template <size_t Order, size_t Layer>
      Matrix_Ref
      compute_neighbour_contribution(const double distance,
                                     const ClusterRefKey<Order, Layer> &) {
        // TODO(alex) TODO(felix) include an check that the distance is within
        // the (x1,x2) range of the interpolator
        this->radial_integral_neighbour = this->intp->interpolate(distance);
        return Matrix_Ref(this->radial_integral_neighbour);
      }

     protected:
      void precompute() override {
        this->precompute_fac_a();
      }
      void precompute_fac_a() {
        auto smearing{downcast_atomic_smearing<AtomicSmearingType::Constant>(
            this->atomic_smearing)};
        this->fac_a = 0.5 * pow(smearing->get_gaussian_sigma(), -2);
      }

      // Should be invoked only after the a-factor has been precomputed
      void init_interpolator(const Hypers_t & hypers) {
        auto radial_contribution_hypers =
            hypers.at("radial_contribution").template get<json>();
        auto optimization_hypers =
            radial_contribution_hypers.at("optimization").template get<json>();

        double accuracy{this->get_interpolator_accuracy(optimization_hypers)};
        double range_begin{this->get_range_begin(optimization_hypers)};
        double range_end{this->get_range_end(optimization_hypers)};
        this->init_interpolator(range_begin, range_end, accuracy);
      }

      void init_interpolator(const double range_begin, const double range_end,
                             const double accuracy) {
        // "this" is passed by reference and is mutable
        std::function<Matrix_t(double)> func{
            [&](const double distance) mutable {
              Parent::compute_neighbour_contribution(distance, this->fac_a);
              return this->radial_integral_neighbour;
            }};
        Matrix_t result = func(range_begin);
        int cols{static_cast<int>(result.cols())};
        int rows{static_cast<int>(result.rows())};
        this->intp = std::make_unique<Interpolator_t>(
            func, range_begin, range_end, accuracy, cols, rows);
      }

      double get_range_begin(const Hypers_t & optimization_hypers) {
        if (optimization_hypers.find("range") != optimization_hypers.end()) {
          return optimization_hypers.at("range")
              .at("begin")
              .template get<double>();
        }
        // default range begin
        return 0.;
      }

      double get_range_end(const Hypers_t & optimization_hypers) {
        if (optimization_hypers.find("range") != optimization_hypers.end()) {
          return optimization_hypers.at("range")
              .at("end")
              .template get<double>();
        }
        throw std::logic_error(
            "Interpolator option is on but no range end for interpolation is "
            "given in the json hyperparameter. Interpolator cannot be "
            "initialized.");
        return 0;
      }

      double get_interpolator_accuracy(const Hypers_t & optimization_hypers) {
        if (optimization_hypers.find("accuracy") != optimization_hypers.end()) {
          return optimization_hypers.at("accuracy").template get<double>();
        }
        // default accuracy
        return 1e-8;
      }

      double get_cutoff(const Hypers_t & hypers) {
        auto fc_hypers = hypers.at("cutoff_function").template get<json>();
        return fc_hypers.at("cutoff").at("value").template get<double>();
      }

      double fac_a{};
      std::unique_ptr<Interpolator_t> intp{};
    };

  }  // namespace internal

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_GEOMETRIC_KSPACE_FACTORS_HH_
