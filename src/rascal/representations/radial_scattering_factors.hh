/**
 * @file   rascal/representations/scattering_factors.hh
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

#ifndef SRC_RASCAL_REPRESENTATIONS_SCATTERING_FACTORS_HH_
#define SRC_RASCAL_REPRESENTATIONS_SCATTERING_FACTORS_HH_

#include "rascal/math/hyp1f1.hh"
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

        // Initialize size of the member arrays
        this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
        this->ortho_norm_matrix.resize(this->max_radial, this->max_radial);
        this->radial_norm_factors.resize(this->max_radial, 1);
        this->radial_sigmas.resize(this->max_radial, 1);
        this->radial_nl_factors.resize(this->max_radial, this->max_angular + 1);
        this->radial_integral.resize(this->max_radial*(this->max_angular + 1));

        // find the cutoff radius of the representation
        auto fc_hypers = hypers.at("cutoff_function").get<json>();
        this->interaction_cutoff = fc_hypers.at("cutoff").at("value").get<double>();

      }

      void precompute() override {
	// compute of radial prefactors
        this->precompute_radial_prefactors();
	// compute orthogonalization matrix
        this->precompute_radial_ortho_matrix();
	// add normalization factors to orthogonalization matrix
        this->ortho_norm_matrix = this->radial_norm_factors.asDiagonal() * this->radial_ortho_matrix;
      }

      /** Function to compute common prefactors for the radial integral */
      void precompute_radial_prefactors() {

        using math::pow;

        for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          // sigmas of radial functions 
          this->radial_sigmas[radial_n] =
              std::max(std::sqrt(static_cast<double>(radial_n)), 1.0) *
              this->interaction_cutoff / static_cast<double>(this->max_radial);
	  // normalization factors
          this->radial_norm_factors(radial_n) =
              std::sqrt(2.0 / (std::tgamma(1.5 + radial_n) *
                        pow(this->radial_sigmas[radial_n], 3 + 2 * radial_n)));
          for (size_t angular_l{0}; angular_l < this->max_angular+1; ++angular_l) {
	     // radial nl-dependent factors  
             this->radial_nl_factors(radial_n,angular_l) =
		 pow(2.0,0.5*(radial_n - angular_l - 1)) 
                 * pow(this->radial_sigmas[radial_n], 3 + radial_n + angular_l)
                 * std::tgamma(0.5 * (3.0 + radial_n + angular_l)) / std::tgamma(1.5+angular_l);
          }
        }
      }

      /** Function to compute the normalized orthogonalization matrix */
      void precompute_radial_ortho_matrix() {

        using math::pow;
        using std::sqrt;
        using std::tgamma;

	// normalized overlap between primitive radial functions
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
                (pow(this->radial_sigmas[radial_n1] 
		     * this->radial_sigmas[radial_n2] , 1.5) *
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
	// orthogonalization matrix
        this->radial_ortho_matrix =
            unitary * eigs_invsqrt.matrix().asDiagonal() * unitary.adjoint();
      }

      // Compute the radial integral I_{nl}(kval)
      Vector_Ref compute_radial_integral(const double kval) {
        
        using math::pow;
	
	double fac_a;
        double fac_b;
        double fac_c;

	int nl_idx{0};
	for (size_t radial_n{0}; radial_n < this->max_radial; ++radial_n) {
          for (size_t angular_l{0}; angular_l < this->max_angular+1; ++angular_l) {

            fac_a = 0.5 * (3.0 + radial_n + angular_l);
	    fac_b = 1.5 + angular_l;
	    fac_c = - 0.5 * pow( kval * this->radial_sigmas[radial_n] , 2); 
            math::Hyp1f1 hyp1f1_calculator(fac_a,fac_b);
            this->radial_integral(nl_idx) = 
		    this->radial_nl_factors(radial_n,angular_l) 
                         * hyp1f1_calculator.calc(fac_c);
	    nl_idx += 1;

          }
        }

        return Vector_Ref(this->radial_integral);

      }

      // Function to normalize and orthogonalize radial projections
      template <typename Coeffs>
      void finalize_coefficients(Coeffs & coefficients) const {
        coefficients.lhs_dot(this->ortho_norm_matrix);
      }

      Hypers_t hypers{};
      // some usefull parameters
      double interaction_cutoff{};
      size_t max_radial{};
      size_t max_angular{};

      Vector_t radial_sigmas{};
      Vector_t radial_norm_factors{};
      Vector_t radial_nl_factors{};
      Matrix_t radial_ortho_matrix{};
      Matrix_t ortho_norm_matrix{};
      Vector_t radial_integral{};

    };

  }  // namespace internal

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_SCATTERING_FACTORS_HH_
