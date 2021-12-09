/**
 * @file   rascal/math/interpolator.hh
 *
 * @author Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 September 2019
 *
 * @brief Implementation of interpolator for functions functions of the form
 *        f:[x1,x2]->ℝ and f:[x1,x2]->ℝ^{n,m}
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
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

#include "rascal/math/interpolator.hh"

using namespace rascal::math;  // NOLINT

bool rascal::math::is_grid_uniform(const Vector_Ref & grid) {
  // checks if the grid is in ascending order
  for (int i = 0; i < grid.size() - 2; i++) {
    if (grid(i + 1) > grid(i)) {
      return false;
    }
  }

  // checks if the cell/step size is the same everywhere
  const double step_size = grid(1) - grid(0);
  for (int i = 0; i < grid.size() - 2; i++) {
    // h_i - h_0 > ε
    if (std::abs((grid(i + 1) - grid(i)) - step_size) > DBL_FTOL) {
      return false;
    }
  }
  return true;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

using CubicSplineScalarUniformInterpolation =
    InterpolationMethod<InterpolationMethod_t::CubicSplineScalarUniform>;

void CubicSplineScalarUniformInterpolation::compute_second_derivative(
    const Vector_Ref & yv, double dfx1, double dfx2) {
  // Bad xa input to routine splint
  int n{static_cast<int>(yv.size())};
  Vector_t y2 = Vector_t::Zero(n);
  Vector_t u = Vector_t::Zero(n);
  double p, qn;

  // forward sweeping
  y2(0) = -0.5;
  u(0) = (3.0 / this->h) * ((yv(1) - yv(0)) / this->h - dfx1);
  for (int i{1}; i < n - 1; i++) {
    p = 0.5 * y2(i - 1) + 2.0;
    y2(i) = -0.5 / p;
    u(i) = (yv(i + 1) - 2 * yv(i) + yv(i - 1)) / this->h;
    u(i) = (3.0 * u(i) / this->h - 0.5 * u(i - 1)) / p;
  }

  // back substitution
  qn = 0.5;
  u(n - 1) = (3.0 / this->h) * (dfx2 - (yv(n - 1) - yv(n - 2)) / this->h);
  y2(n - 1) = (u(n - 1) - qn * u(n - 2)) / (qn * y2(n - 2) + 1.0);
  for (int k{n - 2}; k >= 0; k--) {
    y2(k) = y2(k) * y2(k + 1) + u(k);
  }
  this->second_derivative_h_sq_6 = y2 * this->h_sq_6;
}

void CubicSplineScalarUniformInterpolation::compute_second_derivative(
    const Vector_Ref & yv) {
  // Bad xa input to routine splint
  int n{static_cast<int>(yv.size())};
  Vector_t y2 = Vector_t::Zero(n);
  Vector_t u = Vector_t::Zero(n);
  double p;

  // forward sweeping
  y2(0) = 0.0;
  u(0) = 0.0;
  for (int i{1}; i < n - 1; i++) {
    p = 0.5 * y2(i - 1) + 2.0;
    y2(i) = -0.5 / p;
    u(i) = (yv(i + 1) - 2 * yv(i) + yv(i - 1)) / this->h;
    u(i) = (3.0 * u(i) / this->h - 0.5 * u(i - 1)) / p;
  }

  // back substitution
  y2(n - 1) = 0.0;  // (u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
  for (int k{n - 2}; k > 0; k--) {
    y2(k) = y2(k) * y2(k + 1) + u(k);
  }
  y2(0) = 0.0;  // y2(0)*y2(1)+u(0);
  this->second_derivative_h_sq_6 = y2 * this->h_sq_6;
}

double CubicSplineScalarUniformInterpolation::interpolate_for_one_point(
    const Vector_Ref & xx, const Vector_Ref & yy, int j1, double x) const {
  int klo{j1};
  int khi{j1 + 1};
  // percentage of grid cell start to x
  double a{(xx(khi) - x) / this->h};
  // percentage of x to grid cell end, because a+b = 1 we can simplify
  // the computation
  double b{1 - a};
  return a * (yy(klo) + (a * a - 1) * this->second_derivative_h_sq_6(klo)) +
         b * (yy(khi) + (b * b - 1) * this->second_derivative_h_sq_6(khi));
}

double
CubicSplineScalarUniformInterpolation::interpolate_derivative_for_one_point(
    const Vector_Ref & xx, const Vector_Ref & yy, int j1, double x) const {
  int klo{j1};
  int khi{j1 + 1};
  double a{(xx(khi) - x) / this->h};
  // Because we can assume a+b=1, we simplify the calculation of b
  double b{1 - a};
  return (yy(khi) - yy(klo) -
          (3 * a * a - 1) * this->second_derivative_h_sq_6(klo) +
          (3 * b * b - 1) * this->second_derivative_h_sq_6(khi)) /
         this->h;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

using CubicSplineVectorUniformInterpolation =
    InterpolationMethod<InterpolationMethod_t::CubicSplineVectorUniform>;

void CubicSplineVectorUniformInterpolation::compute_second_derivative(
    const Matrix_Ref & yv) {
  int n{static_cast<int>(yv.rows())};
  Matrix_t y2 = Matrix_t::Zero(n, yv.cols());
  Matrix_t u = Matrix_t::Zero(n, yv.cols());
  double sig = 0.5;
  Vector_t p = Vector_t::Zero(n);
  y2.row(0) = Vector_t::Zero(yv.cols());
  u.row(0) = Vector_t::Zero(yv.cols());
  for (int i{1}; i < n - 1; i++) {
    p = sig * y2.row(i - 1).array() + 2.0;
    y2.row(i) = (sig - 1.0) / p.array();
    // noalias prevents storage of temporaries
    u.row(i).noalias() = (3. *
                          (yv.row(i + 1).array() - 2 * yv.row(i).array() +
                           yv.row(i - 1).array()) /
                          (this->h * this->h))
                             .matrix();
    u.row(i).array() -= sig * u.row(i - 1).array();
    u.row(i).array() /= p.array();
  }
  y2.row(n - 1) = Vector_t::Zero(y2.cols());
  for (int k{n - 2}; k > 0; k--) {
    y2.row(k) = y2.row(k).array() * y2.row(k + 1).array() + u.row(k).array();
  }
  y2.row(0) = Vector_t::Zero(y2.cols());
  this->second_derivative_h_sq_6 = y2 * this->h_sq_6;
}

Vector_t CubicSplineVectorUniformInterpolation::interpolate_for_one_point(
    const Vector_Ref & xx, const Matrix_Ref & yy, int j1, double x) const {
  int klo{j1};
  int khi{j1 + 1};
  // Bad xa input to routine splint
  assert(h > DBL_FTOL);
  double a{(xx(khi) - x) / this->h};
  // Because we can assume a+b=1, we simplify the calculation of b
  double b{1 - a};
  return (a * (yy.row(klo).array() +
               (a * a - 1) * this->second_derivative_h_sq_6.row(klo).array()) +
          b * (yy.row(khi).array() +
               (b * b - 1) * this->second_derivative_h_sq_6.row(khi).array()))
      .matrix();
}

Vector_t
CubicSplineVectorUniformInterpolation::interpolate_derivative_for_one_point(
    const Vector_Ref & xx, const Matrix_Ref & yy, int j1, double x) const {
  int klo{j1};
  int khi{j1 + 1};
  // Bad xa input to routine splint
  assert(h > DBL_FTOL);
  double a{(xx(khi) - x) / this->h};
  // Because we can assume a+b=1, we simplify the calculation of b
  double b{1 - a};
  return ((yy.row(khi).array() - yy.row(klo).array() -
           (3 * a * a - 1) * this->second_derivative_h_sq_6.row(klo).array() +
           (3 * b * b - 1) * this->second_derivative_h_sq_6.row(khi).array()) /
          this->h)
      .matrix();
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void UniformSearchMethod::initialize(const Vector_Ref & grid,
                                     size_t nb_support_points) {
  // nb_grid_points/unit
  this->nb_grid_points_per_unit =
      grid.size() / (grid(grid.size() - 1) - grid(0));
  this->x1 = grid(0);
  this->grid_size = grid.size();
  this->nb_support_points = nb_support_points;
  this->search_size = this->grid_size - this->nb_support_points;
}

int UniformSearchMethod::search(double x, const Vector_Ref &) const {
  return std::max(0, std::min(this->search_size,
                              static_cast<int>((x - this->x1) *
                                               this->nb_grid_points_per_unit) -
                                  1));
}
