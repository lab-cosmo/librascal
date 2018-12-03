/**
 * file   math_utils.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief contains the implementation of miscellaneous math functions
 *
 * Copyright © 2018  Felix Musil, Max Veit, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "math_interface.hh"
#include <Eigen/Dense>
#include <cmath>
#include <limits>
// #include <Eigen/Dense> fails -- WTF?!

namespace rascal {
  namespace math {

    /// Pi to more digits than anyone could possibly need
    const double PI = 3.14159265358979323846264338327950288419716939937510;
    const double SQRT_TWO = std::sqrt((double)2.0);
    const double INV_SQRT_TWO = std::sqrt((double)0.5);
    const double SQRT_THREE = std::sqrt((double)3.0);

    /// How small a number must be to be considered effectively zero
    const double dbl_ftol = 100.0 * std::numeric_limits<double>::epsilon();

    Eigen::MatrixXd compute_assoc_legendre_polynom(
        double cos_theta, size_t max_angular);

    Eigen::MatrixXd compute_cos_sin_angle_multiples(
        double cos_phi, double sin_phi, size_t max_m);

    Eigen::MatrixXd compute_spherical_harmonics(
        const Eigen::Ref<const Eigen::Vector3d> &direction,
        size_t max_angular);
  } // math
} // rascal

#endif /* MATH_UTILS_H */
