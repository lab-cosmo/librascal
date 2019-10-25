/**
 * @file   spherical_harmonics.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 * @author  Alex
 *
 * @date   14 October 2018
 *
 * @brief implementation of the spherical harmonics, optimized, with gradients
 *
 * Copyright  2018  Felix Musil, Max Veit, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_MATH_SPHERICAL_HARMONICS_HH_
#define SRC_MATH_SPHERICAL_HARMONICS_HH_

#include "math_utils.hh"

#include <iostream>
#include <vector>

namespace rascal {
  namespace math {
    /**
     * Compute a full set of spherical harmonics (optimized version)
     *
     * Follows the algorithm described in https://arxiv.org/abs/1410.1748
     *
     * In brief, this class computes the real spherical harmonics including
     * the Condon-Shortley phase where the imaginary components of the usual
     * complex functions are instead stored in the negative-m indices:
     *
     * \f{equation}{
     * Y_\ell^m(\theta, \phi) =
     * \begin{cases}
     * \sqrt{\frac{2\ell + 1}{2\pi} \frac{(\ell + m)!}{(\ell - m)!}}
     *      P_\ell^{-m}(\cos\theta) \sin(-m\phi) & \text{for } m < 0\\
     * \sqrt{\frac{2\ell + 1}{4\pi}} P_\ell(\cos\theta) & \text{for } m = 0\\
     * \sqrt{\frac{2\ell + 1}{2\pi} \frac{(\ell - m)!}{(\ell + m)!}}
     *      P_\ell^{m}(\cos\theta) \cos(m\phi) & \text{for } m > 0
     * \end{cases}
     * \f}
     *
     * In case you're wondering why it's \f$\frac{1}{2π}\f$ on the \f$m \neq
     * 0\f$ components (instead of \f$\frac{1}{4π}\f$), there's an extra factor
     * of 1/2 that comes from integrating cos² or sin² over the full circle of
     * φ, so these are indeed normalized in the same sense as the complex
     * spherical harmonics:
     * \f[
     * \iint_{\mathcal{S}^2} dΩ\; Y_\ell^m(θ, \phi) Y_{\ell'}^{m'}(θ, \phi)
     *      = δ_{\ell\ell'}δ_{mm'}
     * \f]
     *
     * (this extra factor of \f$\frac{1}{\sqrt{2}}\f$ is not included in the
     * normalization of the associated Legendre polynomials defined above;
     * however, all other normalization factors are.)
     *
     * Cartesian gradients can optionally be computed in addition.
     *
     * Part of the efficiency derives from moving the direction-independent
     * calculations into a separate precompute() method; you must therefore call
     * precompute() (which also sets \f$\ell_\text{max}\f$) before any call to
     * calc().
     *
     * Once calc() has been called for a direction vector, the results can be
     * retrieved with get_harmonics() (and the gradients, if computed, with
     * get_harmonics_derivatives())
     *
     * Note the storage order of the harmonics components here follows the
     * "compressed" format, with the \f$\ell\f$ index varying the slowest, from
     * 0 to \f$\ell_\text{max}\f$, and \f$m\f$ varying from \f$-\ell\f$ to
     * \f$\ell\f$ for each value of \f$\ell\f$, for a total of
     * \f$(\ell_\text{max}+1)^2\f$ components.  For example, the first few
     * entries of the array would have the \f$(\ell, m)\f$ numbers: (0,0) (1,-1)
     * (1,0) (1, 1) (2, -2)....
     */
    class SphericalHarmonics {
     protected:
      size_t max_angular{0};
      Vector_t angular_coeffs1{};
      Vector_t angular_coeffs2{};
      Vector_t harmonics{};
      Matrix_t assoc_legendre_polynom{};
      MatrixX2_t cos_sin_m_phi{};
      Matrix_t coeff_a{};
      Matrix_t coeff_b{};
      // derivative related member variables
      bool calculate_derivatives{false};
      bool derivatives_precomputed{false};
      Matrix_t harmonics_derivatives{};
      Matrix_t plm_factors{};
      Vector_t legendre_polynom_differences{};
      Vector_t phi_derivative_factors{};

     public:
      /** Constructs the class, but doesn't precompute or set anything */
      SphericalHarmonics() {}

      /**
       * Construct a SphericalHarmonics class with a default setting for whether
       * to calculate the gradients
       *
       * Remember to call precompute to finish initialization
       */
      explicit SphericalHarmonics(bool calculate_derivatives)
          : calculate_derivatives{calculate_derivatives} {}

      /**
       * Precomputes as many parameters as possible
       *
       * @param max_angular Maximum angular momentum number
       *                    \f$\ell_\text{max}\f$
       *
       * @param calculate_derivatives Whether to precompute for the derivatives
       *     as well.  If the internal default (set at construction) is false,
       *     but the value provided here is true, the internal default will be
       *     set to true.  If the value provided here is false (default
       *     argument) then the internal default is not changed; derivatives
       *     will be precomputed only if the internal default is true.
       */
      void precompute(size_t max_angular, bool calculate_derivatives = false) {
        using math::pow;
        using std::sqrt;
        this->max_angular = max_angular;
        this->harmonics =
            Vector_t::Zero((this->max_angular + 1) * (this->max_angular + 1));
        // Note the ALPs need an extra zero column to allow natural use of the
        // raising and lowering operators (specifically P_l^(l+1), which should
        // give zero).
        this->assoc_legendre_polynom =
            Matrix_t::Zero(this->max_angular + 1, this->max_angular + 2);
        this->coeff_a =
            Matrix_t::Zero(this->max_angular + 1, 2 * this->max_angular + 1);
        this->coeff_b =
            Matrix_t::Zero(this->max_angular + 1, 2 * this->max_angular + 1);
        this->cos_sin_m_phi = MatrixX2_t::Zero(this->max_angular + 1, 2);

        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             angular_l++) {
          double lsq = angular_l * angular_l;
          double lm1sq = (angular_l - 1) * (angular_l - 1);
          // TODO(alex) vectorize
          for (size_t m_count{0}; m_count < angular_l + 1; m_count++) {
            double msq = m_count * m_count;
            this->coeff_a(angular_l, m_count) =
                sqrt((4 * lsq - 1.0) / (lsq - msq));
            this->coeff_b(angular_l, m_count) =
                -1.0 * sqrt((lm1sq - msq) / (4 * lm1sq - 1.0));
          }
        }
        // this->angular_coeffs1.reserve(this->max_angular);
        this->angular_coeffs1 = Vector_t::Zero(this->max_angular + 1);
        this->angular_coeffs2 = Vector_t::Zero(this->max_angular + 1);
        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             angular_l++) {
          angular_coeffs1(angular_l) = sqrt(2 * angular_l + 1);
          angular_coeffs2(angular_l) = -sqrt(1.0 + 0.5 / angular_l);
        }

        // We want to precompute derivative information in almost any case
        if (calculate_derivatives or this->calculate_derivatives) {
          if (not this->calculate_derivatives) {
            this->calculate_derivatives = true;
          }
          this->harmonics_derivatives =
              Matrix_t::Zero(3, pow(this->max_angular + 1, 2));
          this->plm_factors =
              Matrix_t::Zero(this->max_angular + 1, this->max_angular + 1);
          // TODO(alex) can be done with broadcasting and setting half of matrix
          // to zero, but we are actually only calculating the rectangular part,
          // what is faster?
          for (size_t angular_l{0}; angular_l < this->max_angular + 1;
               angular_l++) {
            Vector_t m_counts =
                Eigen::VectorXd::LinSpaced(angular_l + 1, 0, angular_l);
            this->plm_factors.row(angular_l).head(angular_l + 1) =
                ((angular_l - m_counts.array()) *
                 (angular_l + m_counts.array() + 1))
                    .sqrt();
          }
          this->legendre_polynom_differences =
              Vector_t::Zero(this->max_angular);
          this->phi_derivative_factors = Vector_t::Zero(this->max_angular);
          this->derivatives_precomputed = true;
        }
      }

      /**
       * Compute a set of normalized associated Legendre polynomials
       *
       * These are normalized for use in computing the real spherical harmonics;
       * see the class documentation for details.  In particular, the \f$m=0\f$
       * harmonics require an extra factor of \f$\frac{1}{\sqrt{2}}\f$.  The
       * negative-m functions are not computed due to symmetry.
       *
       * @param cos_theta (aka x) Where to evaluate the polynomial
       *
       * Stores the results as an (Eigen)matrix containing the evaluated
       * polynomials, sized \f$\ell_\text{max}\f$ by \f$(2\ell_\text{max} +
       * 1)\f$; the row is indexed by \f$\ell\f$ and the column by \f$m \geq
       * 0\f$.
       */
      void compute_assoc_legendre_polynom(double cos_theta) {
        using math::pow;
        using std::sqrt;
        double sin_theta = sqrt(1.0 - pow(cos_theta, 2));
        const double SQRT_INV_2PI = sqrt(0.5 / PI);
        // Compute the associated Legendre polynomials: l < 2 are special cases
        // These include the normalization factors usually needed in the
        // spherical harmonics
        double l_accum{SQRT_INV_2PI};
        this->assoc_legendre_polynom(0, 0) = SQRT_INV_2PI;
        if (this->max_angular > 0) {
          this->assoc_legendre_polynom(1, 0) =
              cos_theta * SQRT_THREE * SQRT_INV_2PI;
          l_accum = l_accum * -sqrt(3.0 / 2.0) * sin_theta;
          this->assoc_legendre_polynom(1, 1) = l_accum;
        }
        for (size_t angular_l{2}; angular_l < this->max_angular + 1;
             angular_l++) {
          // for l > 1 : Use the recurrence relation
          // TODO(max-veit) don't bother calculating m =/= 0 if sin(theta) == 0
          //                (z-axis)
          // avoid making temp by breaking down the operation in 3 parts
          this->assoc_legendre_polynom.row(angular_l)
              .head(angular_l - 1)
              .array() =
              cos_theta * this->assoc_legendre_polynom.row(angular_l - 1)
                              .head(angular_l - 1)
                              .array();
          this->assoc_legendre_polynom.row(angular_l)
              .head(angular_l - 1)
              .array() +=
              this->coeff_b.row(angular_l).head(angular_l - 1).array() *
              this->assoc_legendre_polynom.row(angular_l - 2)
                  .head(angular_l - 1)
                  .array();
          this->assoc_legendre_polynom.row(angular_l)
              .head(angular_l - 1)
              .array() *=
              this->coeff_a.row(angular_l).head(angular_l - 1).array();

          this->assoc_legendre_polynom(angular_l, angular_l - 1) =
              // cos_theta * sqrt(2 * angular_l + 1) * l_accum;
              l_accum * cos_theta * this->angular_coeffs1(angular_l);
          l_accum = l_accum * sin_theta * this->angular_coeffs2(angular_l);
          this->assoc_legendre_polynom(angular_l, angular_l) = l_accum;
        }
      }

      /**
       * Compute a full set of spherical harmonics given a direction vector.
       * If calculate_derivatives flag is on, the derivatives are additionally
       * computed. This function returns void, results have to be retrieved
       * with get functions.
       *
       * @param direction   unit vector defining the angles
       *                    (arguments for the \f$Y_\ell^m\f$)
       *
       * @param calculate_derivatives       Compute the gradients too?
       *
       * @warning Prints warning and normalizes direction if it is not
       *          already normalized.
       */
      void calc(const Eigen::Ref<const Eigen::Vector3d> & direction,
                bool calculate_derivatives) {
        using math::pow;
        using std::sqrt;
        Eigen::Vector3d direction_normed;
        if (std::abs((direction[0] * direction[0] +
                      direction[1] * direction[1] +
                      direction[2] * direction[2]) -
                     1.0) > math::dbl_ftol) {
          std::cerr << "Warning: SphericalHarmonics::calc()";
          std::cerr << ": Direction vector unnormalized, normalizing it now";
          std::cerr << std::endl;
          direction_normed = direction / direction.norm();
        } else {
          direction_normed = direction;
        }

        // The cosine against the z-axis is just the z-component of the
        // direction vector
        double cos_theta = direction_normed[2];
        // The less efficient, but more intuitive implementation:
        // double phi = std::atan2(direction[1], direction[0]);
        double sqrt_xy = std::hypot(direction_normed[0], direction_normed[1]);
        // For a vector along the z-axis, define phi=0
        double cos_phi{1.0}, sin_phi{0.0};
        if (sqrt_xy >= math::dbl_ftol) {
          cos_phi = direction_normed[0] / sqrt_xy;
          sin_phi = direction_normed[1] / sqrt_xy;
        }

        this->compute_assoc_legendre_polynom(cos_theta);
        this->compute_cos_sin_angle_multiples(cos_phi, sin_phi);

        this->compute_spherical_harmonics();
        if (calculate_derivatives) {
          if (this->derivatives_precomputed) {
            // A rose, by any other name, would have the same exact value
            // double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
            double sin_theta{sqrt_xy};
            this->compute_spherical_harmonics_derivatives(sin_theta, cos_theta,
                                                          sin_phi, cos_phi);
          } else {
            // TODO(max) do we just precompute here instead of throwing a rude
            // error?
            std::stringstream err_str{};
            err_str
                << "Resources for computation of dervatives have not been "
                   "initialized. Please set calculate_derivatives flag on "
                   "construction of the SphericalHarmonics object or during "
                   "precomputation.";
            throw std::runtime_error(err_str.str());
          }
        }
      }

      /**
       * Same as calc(), but using the internal default to decide whether to
       * compute derivatives.
       */
      void calc(const Eigen::Ref<const Eigen::Vector3d> & direction) {
        this->calc(direction, this->calculate_derivatives);
      }

      /**
       * Compute \f$\cos(m\phi)\f$ and \f$\sin(m\phi)\f$ from the recurrence
       * relations
       *
       * The relations are (these are the same recurrence relations used to
       * calculate the Chebyshev polynomials):
       *
       * \f{eqnarray}{
       * \cos(m\phi) &=& 2\cos(\phi)\cos((m-1)\phi) - \cos((m-2)\phi)\\
       * \sin(m\phi) &=& 2\cos(\phi)\sin((m-1)\phi) - \sin((m-2)\phi)
       * \f}
       *
       * and they require only \f$\cos(\phi)\f$ and \f$\sin(\phi)\f$ to start.
       *
       * @param cos_phi Value of \f$\cos(\phi)\f$ to start the relation
       *
       * @param sin_phi Value of \f$\sin(\phi)\f$ to start the relation
       *
       * Stores the results as an (Eigen)matrix, sized \f$m_\text{max}\f$ by 2
       * with the \f$\cos(m\phi)\f$ stored in the first column and
       * \f$\sin(m\phi)\f$ in the second column, with \f$m\f$ being the row
       * index
       */
      void compute_cos_sin_angle_multiples(double cos_phi, double sin_phi) {
        for (size_t m_count{0}; m_count < this->max_angular + 1; m_count++) {
          if (m_count == 0) {
            this->cos_sin_m_phi.row(m_count) << 1.0, 0.0;
          } else if (m_count == 1) {
            this->cos_sin_m_phi.row(m_count) << cos_phi, sin_phi;
          } else {
            this->cos_sin_m_phi.row(m_count) =
                2.0 * cos_phi * this->cos_sin_m_phi.row(m_count - 1) -
                this->cos_sin_m_phi.row(m_count - 2);
          }
        }
      }

      /**
       * Combine the associated Legendre polynomials with the cos/sin(mφ) to
       * compute the final spherical harmonics
       *
       * Stores the results as an (Eigen)matrix, sized
       * \f$(\ell_\text{max}+1)^2\f$.  The index collects \f$\ell\f$ and \f$m\f$
       * quantum numbers stored in compact format (\f$m\f$ varies fastest, from
       * \f$-\ell\f$ to \f$\ell\f$, and \f$\ell\f$ from 0 to
       * \f$\ell_\text{max}\f$).
       */
      void compute_spherical_harmonics() {
        size_t lm_base{0};  // starting point for storage
        for (size_t angular_l{0}; angular_l < this->max_angular + 1;
             angular_l++) {
          // uses symmetry of spherical harmonics,
          // careful with the storage order
          // TODO(alex) please clarify -- is it the same storage order we had
          // before (i.e. negative-m components first)?  If so, please refer to
          // the above documentation or delete this comment -- Max
          this->harmonics.segment(lm_base + angular_l, angular_l + 1) =
              this->assoc_legendre_polynom.row(angular_l)
                  .head(angular_l + 1)
                  .array() *
              cos_sin_m_phi.col(0).head(angular_l + 1).transpose().array();
          this->harmonics.segment(lm_base, angular_l + 1) =
              (this->assoc_legendre_polynom.row(angular_l)
                   .head(angular_l + 1)
                   .array() *
               cos_sin_m_phi.col(1).head(angular_l + 1).transpose().array())
                  .reverse();
          this->harmonics(lm_base + angular_l) =
              this->assoc_legendre_polynom(angular_l, 0) * INV_SQRT_TWO;
          lm_base += 2 * angular_l + 1;
        }  // for (l in [0, lmax])
      }

      /**
       * Compute a full set of spherical harmonics and their Cartesian gradients
       *
       * Spherical harmonics are defined as described in
       * compute_spherical_harmonics().  Gradients are defined with respect
       * to motion of the central atom, which is the opposite sign of the usual
       * definition with respect to the _arguments_ of the \f$Y_\ell^m\f$.  The
       * actual Cartesian gradients include an extra factor of \f$\frac{1}{r}\f$
       * that is not included here; the rest is independent of radius.
       *
       * Results are stored as in compute_spherical_harmonics(), but with an
       * extra (row) index for the x, y, and z Cartesian components.
       */
      void compute_spherical_harmonics_derivatives(double sin_theta,
                                                   double cos_theta,
                                                   double sin_phi,
                                                   double cos_phi) {
        using std::pow;
        using std::sqrt;

        // angular_l = 0
        // d/dx, d/dy, d/dz
        this->harmonics_derivatives.col(0).setZero();

        // angular_l > 0
        size_t l_block_index{1};
        for (size_t angular_l{1}; angular_l < this->max_angular + 1;
             angular_l++) {
          // legendre_polynom_difference
          legendre_polynom_differences.head(angular_l) =
              this->plm_factors.row(angular_l).segment(0, angular_l).array() *
                  this->assoc_legendre_polynom.row(angular_l)
                      .segment(0, angular_l)
                      .array() -
              this->plm_factors.row(angular_l).segment(1, angular_l).array() *
                  this->assoc_legendre_polynom.row(angular_l)
                      .segment(2, angular_l)
                      .array();
          // TODO(alex) this if then else could be optimized, but I guess the
          // compiler does it.
          // phi_derivative_factor
          if (sin_theta > 0.1) {
            // TODO(alex) could be precomputed
            Vector_t m_counts =
                Eigen::VectorXd::LinSpaced(angular_l, 1, angular_l);
            // singularity at the poles
            phi_derivative_factors.head(angular_l) =
                m_counts.head(angular_l).array() / sin_theta *
                this->assoc_legendre_polynom.row(angular_l)
                    .segment(1, angular_l)
                    .array();
          } else {
            // singularity at the equator
            phi_derivative_factors.head(angular_l) =
                -0.5 / cos_theta *
                (this->plm_factors.row(angular_l)
                         .segment(0, angular_l)
                         .array() *
                     this->assoc_legendre_polynom.row(angular_l)
                         .segment(0, angular_l)
                         .array() +
                 this->plm_factors.row(angular_l)
                         .segment(1, angular_l)
                         .array() *
                     this->assoc_legendre_polynom.row(angular_l)
                         .segment(2, angular_l)
                         .array());
          }

          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //                            ‾‾‾‾‾‾‾

          // d/dx
          this->harmonics_derivatives(0, l_block_index + angular_l) =
              cos_theta * cos_phi * this->plm_factors(angular_l, 0) *
              INV_SQRT_TWO * this->assoc_legendre_polynom(angular_l, 1);
          // d/dy
          this->harmonics_derivatives(1, l_block_index + angular_l) =
              cos_theta * sin_phi * this->plm_factors(angular_l, 0) *
              INV_SQRT_TWO * this->assoc_legendre_polynom(angular_l, 1);
          // d/dz
          this->harmonics_derivatives(2, l_block_index + angular_l) =
              -1.0 * sin_theta * this->plm_factors(angular_l, 0) *
              INV_SQRT_TWO * this->assoc_legendre_polynom(angular_l, 1);

          // d/dx
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //                                    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(0).segment(
              l_block_index + angular_l + 1, angular_l) =
              sin_phi * phi_derivative_factors.head(angular_l).array() *
              this->cos_sin_m_phi.col(1)
                  .segment(1, angular_l)
                  .transpose()
                  .array();

          this->harmonics_derivatives.row(0)
              .segment(l_block_index + angular_l + 1, angular_l)
              .array() += -0.5 * cos_theta * cos_phi *
                          this->cos_sin_m_phi.col(0)
                              .segment(1, angular_l)
                              .transpose()
                              .array() *
                          legendre_polynom_differences.head(angular_l).array();
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(0)
              .segment(l_block_index, angular_l)
              .reverse() = -1.0 * sin_phi *
                           phi_derivative_factors.head(angular_l).array() *
                           this->cos_sin_m_phi.col(0)
                               .segment(1, angular_l)
                               .transpose()
                               .array();
          this->harmonics_derivatives.row(0)
              .segment(l_block_index, angular_l)
              .reverse()
              .array() += -0.5 * cos_theta * cos_phi *
                          this->cos_sin_m_phi.col(1)
                              .segment(1, angular_l)
                              .transpose()
                              .array() *
                          legendre_polynom_differences.head(angular_l).array();

          // d/dy
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //                                    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(1).segment(
              l_block_index + angular_l + 1, angular_l) =
              -1.0 * cos_phi * phi_derivative_factors.head(angular_l).array() *
              this->cos_sin_m_phi.col(1)
                  .segment(1, angular_l)
                  .transpose()
                  .array();
          this->harmonics_derivatives.row(1)
              .segment(l_block_index + angular_l + 1, angular_l)
              .array() += -0.5 * cos_theta * sin_phi *
                          this->cos_sin_m_phi.col(0)
                              .segment(1, angular_l)
                              .transpose()
                              .array() *
                          legendre_polynom_differences.head(angular_l).array();
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(1)
              .segment(l_block_index, angular_l)
              .reverse() = cos_phi *
                           phi_derivative_factors.head(angular_l).array() *
                           this->cos_sin_m_phi.col(0)
                               .segment(1, angular_l)
                               .transpose()
                               .array();
          this->harmonics_derivatives.row(1)
              .segment(l_block_index, angular_l)
              .reverse()
              .array() += -0.5 * cos_theta * sin_phi *
                          this->cos_sin_m_phi.col(1)
                              .segment(1, angular_l)
                              .transpose()
                              .array() *
                          legendre_polynom_differences.head(angular_l).array();
          // d/dz
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //                                    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(2).segment(
              l_block_index + angular_l + 1, angular_l) =
              0.5 * sin_theta *
              this->cos_sin_m_phi.col(0)
                  .segment(1, angular_l)
                  .transpose()
                  .array() *
              legendre_polynom_differences.head(angular_l).array();
          // (l_block - angular_l,... , l_block, ..., l_block + angular_l)
          //  ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          this->harmonics_derivatives.row(2)
              .segment(l_block_index, angular_l)
              .reverse() = 0.5 * sin_theta *
                           this->cos_sin_m_phi.col(1)
                               .segment(1, angular_l)
                               .transpose()
                               .array() *
                           legendre_polynom_differences.head(angular_l).array();

          l_block_index += (2 * angular_l + 1);
        }  // for (l in [0, lmax])
      }

      const MatrixX2_Ref get_cos_sin_m_phi() {
        return MatrixX2_Ref(this->cos_sin_m_phi);
      }

      // Since for calculation purposes assoc_legendre_polynom has one column
      // more than it would have in standard libaries, we return only the
      // segment of size (max_angular+1, max_angular+1) as other standard
      // libaries.
      const Matrix_Ref get_assoc_legendre_polynom() {
        return Matrix_Ref(this->assoc_legendre_polynom)
            .topLeftCorner(this->max_angular + 1, this->max_angular + 1);
      }

      // Returns the (max_angular+1, max_angular+2) matrix related to the
      // associated legendre polynomial with only zeros entries in the last
      // column. Used for testing purposes.
      const Matrix_Ref get_assoc_legendre_polynom_raw() {
        return Matrix_Ref(this->assoc_legendre_polynom);
      }

      /**
       * Access the computed spherical harmonics.
       *
       * @return  (Eigen)matrix containing the results.
       *          Sized \f$(\ell_\text{max}+1)^2\f$, the index collects
       *          \f$\ell\f$ and \f$m\f$ quantum numbers stored in compact
       *          format (\f$m\f$ varies fastest, from \f$-\ell\f$ to
       *          \f$\ell\f$, and \f$\ell\f$ from 0 to \f$\ell_\text{max}\f$).
       */
      const Vector_Ref get_harmonics() { return Vector_Ref(this->harmonics); }

      /**
       * Access the Cartesian gradients of the harmonics, if computed
       *
       * @todo (alex) what happens if they haven't been computed and this gets
       *              called anyway?
       *
       * @return  (Eigen)matrix containing the results.
       *          Sized 3 by \f$(\ell_\text{max}+1)^2\f$; the first index runs
       *          over the the x, y, and z gradient components and the second
       *          index collects \f$\ell\f$ and \f$m\f$ quantum numbers, stored
       *          in compact format (\f$m\f$ varies fastest, from \f$-\ell\f$ to
       *          \f$\ell\f$, and \f$\ell\f$ from 0 to \f$\ell_\text{max}\f$).
       */
      const Matrix_Ref get_harmonics_derivatives() {
        return Matrix_Ref(this->harmonics_derivatives);
      }
    };

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_SPHERICAL_HARMONICS_HH_
