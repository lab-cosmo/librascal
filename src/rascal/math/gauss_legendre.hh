/**
 * @file   rascal/math/gauss_legendre.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   27 May 2019
 *
 * @brief Implementation of the gauss legendre quadrature weights and points
 *
 * Copyright  2019  Felix Musil, Max Veit COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef SRC_RASCAL_MATH_GAUSS_LEGENDRE_HH_
#define SRC_RASCAL_MATH_GAUSS_LEGENDRE_HH_

#include "rascal/math/utils.hh"

namespace rascal {
  namespace math {
    /**
     * Produce the weights and points for the Gauss-Legendre quadrature.
     *
     * @param r_st starting point of the integral
     * @param r_nd ending point of the integral
     * @param order_n number of gauss-legendre quadrature points
     *
     * @return point_weight order_n x 2 matrix containing the points (1st
     * column) and the weights (2nd column) of the Gauss-Legendre quadrature.
     *
     * see https://github.com/JuliaApproximation/FastGaussQuadrature.jl
     * for potential improvements of the method used here (Numerical Recepies)
     */
    MatrixX2_t compute_gauss_legendre_points_weights(double r_st, double r_nd,
                                                     int order_n);

  }  // namespace math
}  // namespace rascal

#endif  // SRC_RASCAL_MATH_GAUSS_LEGENDRE_HH_
