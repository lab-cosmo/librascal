/**
 * file   spherical_harmonics.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief implementation of the spherical harmonics
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
                                            size_t max_angular);

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
                                               size_t max_m);

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
     * @return  (Eigen)vector containing the results.
     *          Sized (l_max+1)**2, contains the l,m components in compressed
     *          format, i.e. (00)(1-1)(10)(11)(2-2)....
     */
    Vector_t compute_spherical_harmonics(
        const Eigen::Ref<const Eigen::Vector3d> & direction,
        size_t max_angular);

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_SPHERICAL_HARMONICS_HH_
