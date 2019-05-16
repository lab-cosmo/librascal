/**
 * file   math_utils.cc
 *

 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief contains the implementation of miscellaneous math functions
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

#include "math_utils.hh"

namespace rascal {
  namespace math {

    /**
     * Compute a set of normalized associated Legendre polynomials
     *
     * These are normalized for use in computing the real spherical harmonics;
     * see math::compute_spherical_harmonics() for details.  The m==0 harmonics
     * require an extra factor of 1/√2 (see below).  Only positive-m functions
     * are computed.
     *
     * @param cos_theta (aka x) Where to evaluate the polynomial
     *
     * @param max_angular Compute up to this angular momentum number (l_max)
     *
     * @return assoc_legendre_polynom
     *        An (Eigen)matrix containing the evaluated polynomials.
     *        Sized l_max by (2*lmax + 1); the row is indexed by l and the
     *        column by m >= 0.
     */
    Matrix_t compute_assoc_legendre_polynom(double cos_theta,
                                            size_t max_angular) {
      using math::pow;
      using std::sqrt;
      /// Technically abs(sin(θ)), but θ only goes from [0, π)
      double sin_theta = sqrt(1.0 - pow(cos_theta, 2));
      // Matrix_t assoc_legendre_polynom(max_angular + 1, max_angular +
      // 1);
      Matrix_t assoc_legendre_polynom =
          Matrix_t::Zero(max_angular + 1, max_angular + 1);
      Matrix_t coeff_a = Matrix_t::Zero(max_angular + 1, 2 * max_angular + 1);
      Matrix_t coeff_b = Matrix_t::Zero(max_angular + 1, 2 * max_angular + 1);
      const double SQRT_INV_2PI = sqrt(0.5 / PI);

      // Coefficients for computing the associated Legendre polynomials
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        double lsq = angular_l * angular_l;
        double lm1sq = (angular_l - 1) * (angular_l - 1);
        for (size_t m_count{0}; m_count < angular_l + 1; m_count++) {
          double msq = m_count * m_count;
          coeff_a(angular_l, m_count) = sqrt((4 * lsq - 1.0) / (lsq - msq));
          coeff_b(angular_l, m_count) =
              -1.0 * sqrt((lm1sq - msq) / (4 * lm1sq - 1.0));
        }
      }
      // Compute the associated Legendre polynomials: l < 2 are special cases
      // These include the normalization factors usually needed in the spherical
      // harmonics
      double l_accum{SQRT_INV_2PI};
      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        if (angular_l == 0) {
          assoc_legendre_polynom(angular_l, 0) = SQRT_INV_2PI;
          continue;
        } else if (angular_l == 1) {
          assoc_legendre_polynom(angular_l, 0) =
              cos_theta * SQRT_THREE * SQRT_INV_2PI;
          l_accum = l_accum * -1.0 * sqrt(3.0 / 2.0) * sin_theta;
          assoc_legendre_polynom(angular_l, 1) = l_accum;
          continue;
        }
        // for l > 1 : Use the recurrence relation
        // TODO(max-veit) don't bother calculating m =/= 0 if sin(theta) == 0
        //                (z-axis)
        for (size_t m_count{0}; m_count < angular_l - 1; m_count++) {
          assoc_legendre_polynom(angular_l, m_count) =
              coeff_a(angular_l, m_count) *
              (cos_theta * assoc_legendre_polynom(angular_l - 1, m_count) +
               coeff_b(angular_l, m_count) *
                   assoc_legendre_polynom(angular_l - 2, m_count));
        }
        assoc_legendre_polynom(angular_l, angular_l - 1) =
            cos_theta * sqrt(2.0 * (angular_l - 1) + 3) * l_accum;
        l_accum = l_accum * sin_theta * -1.0 * sqrt(1.0 + 0.5 / angular_l);
        assoc_legendre_polynom(angular_l, angular_l) = l_accum;
      }
      return assoc_legendre_polynom;
    }

    /**
     * Compute cos(mφ) and sin(mφ) from the recurrence relations
     *
     * The relations are (these are the same recurrence relations used to
     * calculate the Chebyshev polynomials):
     *
     * cos(mφ) = 2cos(φ)cos((m-1)φ) - cos((m-2)φ)
     * sin(mφ) = 2cos(φ)sin((m-1)φ) - sin((m-2)φ)
     *
     * and they require only cos(φ) and sin(φ) to start.
     *
     * @param cos_phi Value of cos(φ) to start the relation
     *
     * @param sin_phi Value of sin(φ) to start the relation
     *
     * @param max_m Compute up to a maximum value of max_m (inclusive)
     *
     * @return cos_sin_m_phi
     *        (Eigen)matrix containing the results.
     *        Sized max_m by 2 with the cos(mφ) stored in the first column
     *        and sin(mφ) in the second column, m being the row index
     */
    MatrixX2_t compute_cos_sin_angle_multiples(double cos_phi, double sin_phi,
                                               size_t max_m) {
      MatrixX2_t cos_sin_m_phi = MatrixX2_t::Zero(max_m + 1, 2);
      for (size_t m_count{0}; m_count < max_m + 1; m_count++) {
        if (m_count == 0) {
          cos_sin_m_phi.row(m_count) << 1.0, 0.0;
        } else if (m_count == 1) {
          cos_sin_m_phi.row(m_count) << cos_phi, sin_phi;
        } else {
          cos_sin_m_phi.row(m_count) =
              2.0 * cos_phi * cos_sin_m_phi.row(m_count - 1) -
              cos_sin_m_phi.row(m_count - 2);
        }
      }
      return cos_sin_m_phi;
    }

    /**
     * Compute a full set of spherical harmonics given a direction vector
     *
     * Follows the algorithm described in https://arxiv.org/abs/1410.1748
     *
     * In brief, this computes the real spherical harmonics where the imaginary
     * components of the usual complex functions are instead stored in the
     * negative-m indices:
     *
     *               ╭ √((2l+1)/(2*π) * (l+m)!/(l-m)!) P_l^-m(cos(θ)) sin(-mφ)
     *               |                                                  for m<0
     * Y_l^m(θ, φ) = ┤ √((2l+1)/(4*π)) P_l(cos(θ)) for m==0
     *               |
     *               ╰ √((2l+1)/(2*π) * (l-m)!/(l+m)!) P_l^m(cos(θ)) cos(mφ)
     *                                                                  for m>0
     *
     * In case you're wondering why it's 1/2π on the m=/=0 components (instead
     * of 1/4π), there's an extra factor of 1/2 that comes from integrating cos²
     * or sin² over the full circle of φ, so these are indeed normalized in the
     * same sense as the complex spherical harmonics:
     *
     * ∫∫_(Sphere) dΩ Y_l^m(θ, φ) Y_l'^m'(θ, φ) = δ_(ll')δ_(mm')
     *
     * (this extra factor of 1/√2 is not included in the normalization of the
     * associated Legendre polynomials defined above; however, all other
     * normalization factors are.)
     *
     * @param direction Unit vector giving the angles (arguments for the Y_l^m)
     *
     * @param max_angular Compute up to this angular momentum number (l_max)
     *
     * @return  (Eigen)matrix containing the results.
     *          Sized l_max+1 by 2l_max+1, the row index
     *          corresponding to l numbers and the column index to m.
     */
    Matrix_t compute_spherical_harmonics(
        const Eigen::Ref<const Eigen::Vector3d> & direction,
        size_t max_angular) {
      using math::pow;
      using std::sqrt;

      if (direction.size() != 3) {
        throw std::length_error("Direction must be a vector in R^3");
      }
      // The cosine against the z-axis is just the z-component of the
      // direction vector
      double cos_theta = direction[2];
      // The less efficient, but more intuitive implementation:
      // double phi = std::atan2(direction[1], direction[0]);
      double sqrt_xy = std::hypot(direction[0], direction[1]);
      // For a vector along the z-axis, define phi=0
      double cos_phi{1.0}, sin_phi{0.0};
      if (sqrt_xy >= math::dbl_ftol) {
        cos_phi = direction[0] / sqrt_xy;
        sin_phi = direction[1] / sqrt_xy;
      }
      Matrix_t harmonics = Matrix_t::Zero(max_angular + 1, 2 * max_angular + 1);
      Matrix_t assoc_legendre_polynom =
          compute_assoc_legendre_polynom(cos_theta, max_angular);
      MatrixX2_t cos_sin_m_phi =
          compute_cos_sin_angle_multiples(cos_phi, sin_phi, max_angular);

      for (size_t angular_l{0}; angular_l < max_angular + 1; angular_l++) {
        for (size_t m_count{0}; m_count < angular_l + 1; m_count++) {
          if (m_count == 0) {
            harmonics(angular_l, angular_l) =
                assoc_legendre_polynom(angular_l, m_count) * INV_SQRT_TWO;
          } else {
            harmonics(angular_l, angular_l + m_count) =
                assoc_legendre_polynom(angular_l, m_count) *
                cos_sin_m_phi(m_count, 0);
            harmonics(angular_l, angular_l - m_count) =
                assoc_legendre_polynom(angular_l, m_count) *
                cos_sin_m_phi(m_count, 1);
          }  // if (m_count == 0)
        }    // for (m_count in [0, l])
      }      // for (l in [0, lmax])
      return harmonics;
    }  // compute_spherical_harmonics()

  }  // namespace math
}  // namespace rascal
