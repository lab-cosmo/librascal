/**
 * file   test_sph_harm_poles.cc
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   13 May 2019
 *
 * @brief  Numerical test of spherical harmonics gradient formule
 *
 * Copyright © 2018 Max Veit, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "math/math_utils.hh"
#include "math/spherical_harmonics.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

using rascal::math::PI;
using rascal::math::SphericalHarmonics;

int main() {
  double cos_theta{0.};
  double sin_theta{0.};
  size_t m_count{1};
  size_t angular_l{2};
  double lowering_plm_factor{
      sqrt((angular_l + m_count) * (angular_l - m_count + 1))};
  double raising_plm_factor{
      sqrt((angular_l - m_count) * (angular_l + m_count + 1))};
  Eigen::MatrixXd assoc_legendre_polynom;
  SphericalHarmonics harmonics_calculator{false};
  harmonics_calculator.precompute(angular_l);
  std::cout << "l = " << angular_l << "\tm = " << m_count << std::endl;
  std::cout << "# theta\tPole singularity\tEquator singularity" << std::endl;
  std::cout << std::setprecision(15);
  for (double theta{PI / 4.}; theta > 1E-15; theta = theta / 2.) {
    cos_theta = std::cos(theta);
    sin_theta = std::sqrt(1.0 - pow(cos_theta, 2));
    harmonics_calculator.compute_assoc_legendre_polynom(cos_theta);
    assoc_legendre_polynom = harmonics_calculator.get_assoc_legendre_polynom();
    double result_nopoles{m_count * assoc_legendre_polynom(angular_l, m_count) /
                          sin_theta};
    double result_noequator{
        -0.5 / cos_theta *
        (lowering_plm_factor * assoc_legendre_polynom(angular_l, m_count - 1) +
         raising_plm_factor * assoc_legendre_polynom(angular_l, m_count + 1))};
    std::cout << theta << "\t";
    std::cout << result_nopoles << "\t";
    std::cout << result_noequator << std::endl;
  }
  for (double theta_c{PI / 4.}; theta_c > 1E-15; theta_c = theta_c / 2.) {
    sin_theta = std::cos(theta_c);
    cos_theta = std::sqrt(1.0 - pow(sin_theta, 2));
    harmonics_calculator.compute_assoc_legendre_polynom(cos_theta);
    assoc_legendre_polynom = harmonics_calculator.get_assoc_legendre_polynom();
    double result_nopoles{m_count * assoc_legendre_polynom(angular_l, m_count) /
                          sin_theta};
    double result_noequator{
        -0.5 / cos_theta *
        (lowering_plm_factor * assoc_legendre_polynom(angular_l, m_count - 1) +
         raising_plm_factor * assoc_legendre_polynom(angular_l, m_count + 1))};
    std::cout << PI / 2. - theta_c << "\t";
    std::cout << result_nopoles << "\t";
    std::cout << result_noequator << std::endl;
  }
}
