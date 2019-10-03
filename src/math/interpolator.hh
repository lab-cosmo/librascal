/**
 * file   interpolator.hh
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

#ifndef SRC_MATH_INTERPOLATOR_HH_
#define SRC_MATH_INTERPOLATOR_HH_

#include <functional>
#include <forward_list>
#include <iostream>
#include <limits>
#include <cassert>
#include "math_utils.hh"

namespace rascal {
  namespace math {
    /**
     * To allow flexibility for the behaviour of the interpolator for
     * experimenting while giving the possibility to optimize the interpolator
     * for certain methods, the interpolator class itself is only a skeleton to
     * execute these interpolator rationals/methods.
     *
     * An interpolator can be constructed given an error bound or a uniform
     * grid. The error is computed by creating a `test_grid` and compare the
     * function results with the results of the interpolator. If the error is
     * not satisfied the grid is refined and the procedure is repeated.
     * 
     * Legend for documentation:
     * f,y         the function to be interpolated
     * Ω           the function is declared as f:[x1,x2]->Ω, where Ω can be ℝ or 
     *             ℝ^{n,m}
     * intp, s     the interpolator function
     * grid        the grid for interpolation
     * test_grid   the test grid to estimate the error
     * ε           error function 
     *
     *  
     * The interpolation method is
     * contained in the InterpolationMethod, for the creation of the grid and
     * test grid contained in the GridRational, a search method to find the
     * nearest grid point before a point for interpolation contained in
     * SearchMethod and a method for computing an error term contained in
     * ErrorMethod.
     *
     * The Interpolator class hierarchy is split into InterpolatorBase ->
     * Interpolator (abstract) -> YOUR_INTERPOLATOR_IMPLEMENTATION An
     * interpolator implementation can be very general an allow multiple
     * combination of methods to allow flexibility or very rigid to allow better
     * optimization. The same goes for the interpolator methods which can be
     * used interchangeable with other methods or depend that certain other
     * methods are used.
     *
     * While testing different methods it has turned out that the interpolator
     * is most optimal for uniform grids, therefore the interpolation methods
     * only work with uniform grids and the only supported grid type is uniform.
     *
     * There has been an implementation for an adaptive refinement method and
     * search methods for non adaptive grids, which were removed but can be
     * found here if needed:
     * https://gist.github.com/agoscinski/59b103998c7faae979dbbc672e8048f6.js
     */

    enum class RefinementMethod_t { Exponential, Linear };

    /**
     * GridRational creates the grid and test grid. The test grid is used
     * to estimate the error of the interpolator on the grid. If the error is
     * too large, the grid is further refined until the desired accuracy is
     * achieved. GridRational consists of an refinement part which determines
     * the method of refining the grid.
     * 
     */
    template <RefinementMethod_t Method>
    struct GridRational {
      /**
       * Computes a grid for interpolation.
       *
       * @param x1 begin of range
       * @param x2 end of range
       * @param degree_of_fineness determines the number of grid points, has to
       *        be greater than 1.
       */
      Vector_t compute_grid(double x1, double x2, int degree_of_fineness) {
        throw std::runtime_error("Function `compute_grid` has not been "
                                 "implemented in GridRational.");
        return 0;
      }

      /**
       * Computes a test grid to estimate the error on the `grid`. The test grid
       * is always the points between the grid points
       *
       * @param x1 begin of range
       * @param x2 end of range
       * @param degree_of_fineness determines the number of grid points, has to
       *        be greater than 1.
       */
      Vector_t compute_test_grid(double x1, double x2, int degree_of_fineness) {
        throw std::runtime_error("Function `compute_test_grid` has not been "
                                 "implemented in GridRational.");
        return 0;
      }
    };

    /**
     * Linear GridRational increases the number of points by `slope` times
     * `degree_of_fineness` in each step.
     *
     * The uniform grid uses the points in between grid points as test grid. The
     * refined grid adds `slope` number of points to grid.
     *
     *                x1                                      x2
     * grid           [       |       |       |       |       ]
     * test grid          x       x       x       x       x
     * refined grid   [    |    |    |    |    |    |    |    ] for slope = 3
     */
    template <>
    class GridRational<RefinementMethod_t::Linear> {
     public:
      explicit GridRational<RefinementMethod_t::Linear>(
          int slope = 10)
          : slope{slope} {}

      Vector_t compute_grid(double x1, double x2, int degree_of_fineness) {
        double nb_grid_points = degree_of_fineness * slope + 2;
        return Vector_t::LinSpaced(nb_grid_points, x1, x2);
      }
      Vector_t compute_test_grid(double x1, double x2, int degree_of_fineness) {
        double nb_grid_points = degree_of_fineness * slope + 2;
        // half of a step size in the grid to get the inner grid points
        // step size = (x2-x1)/(nb_grid_points-1)
        double offset{(x2 - x1) / (2 * (nb_grid_points - 1))};
        return Vector_t::LinSpaced(nb_grid_points - 1, x1 + offset,
                                   x2 - offset);
      }
      bool is_grid_uniform(const Vector_Ref & grid) {
        double step_size = grid(1) - grid(0);
        for (int i = 0; i < grid.size() - 2; i++) {
          // h_i - h_0 > ε
          if (std::abs((grid(i + 1) - grid(i)) - step_size) > dbl_ftol) {
            return false;
          }
        }
        return true;
      }

     private:
      const int slope;
    };

    /**
     * Exponential GridRational increases the number of points by
     * pow(`base`,`degree_of_fineness`) in each step.
     *
     *              x1                  x2
     * grid         [   |   |   |   |   ]
     * test grid      x   x   x   x   x
     * refined grid [ | | | | | | | | | ]  for base = 2
     */
    template <>
    class GridRational<RefinementMethod_t::Exponential> {
     public:
      GridRational<RefinementMethod_t::Exponential>(
          int base=2)
          : base{base} {}

      // TODO(alex) use base
      Vector_t compute_grid(double x1, double x2, int degree_of_fineness) {
        double nb_grid_points = 2 << degree_of_fineness;
        return Vector_t::LinSpaced(nb_grid_points, x1, x2);
      }
      Vector_t compute_test_grid(double x1, double x2, int degree_of_fineness) {
        double nb_grid_points = 2 << degree_of_fineness;
        // half of a step size in the grid to get the inner grid points
        // step size = (x2-x1)/(nb_grid_points-1)
        double offset{(x2 - x1) / (2 * (nb_grid_points - 1))};
        return Vector_t::LinSpaced(nb_grid_points - 1, x1 + offset,
                                   x2 - offset);
      }
      bool is_grid_uniform(const Vector_Ref & grid) {
        double step_size = grid(1) - grid(0);
        for (int i = 0; i < grid.size() - 2; i++) {
          // h_i - h_0 > ε
          if (std::abs((grid(i + 1) - grid(i)) - step_size) > dbl_ftol) {
            return false;
          }
        }
        return true;
      }
      void update_errors(Vector_Ref);

     private:
      const int base;
    };

    enum class InterpolationMethod_t {
      CubicSplineScalarUniform,
      CubicSplineVectorUniform
    };

    /**
     * Interpolation method for a function of the form f:[x1,x2]->Ω
     */
    template <InterpolationMethod_t Type>
    struct InterpolationMethod {
      /**
       * Initializes or reinitializes the interpolation method for a new grid
       * for a functions of the form f:[x1,x2]->ℝ
       *
       * @param grid the grid for interpolation
       * @param evaluated_grid evaluation of the function f(grid)
       */
      void initialize(const Vector_Ref & grid,
                      const Vector_Ref & evaluated_grid) {
        throw std::runtime_error(
            "Function `initialize` for functions of the form f:[x1,x2]->ℝ has "
            "not been implemented in InterpolationMethod.");
       }

      /**
       * Initializes or reinitializes the interpolation method for a new grid
       * for a function of the form f:[x1,x2]->ℝ^{n,m}
       *
       * @param grid the grid for interpolation
       * @param evaluated_grid evaluation of the function f(grid) of shape
       *                       (grid_size, n*m)
       */
      void initialize(const Vector_Ref & grid,
                     const Matrix_Ref & evaluated_grid) {
        throw std::runtime_error(
            "Function `initialize` for functions of the form f:[x1,x2]->ℝ^{n,m} has "
            "not been implemented in InterpolationMethod.");
      }

      /**
       * intp(x) for functions of the form f:[x1,x2]->ℝ.
       *
       * @param grid the grid for interpolation
       * @param evaluated_grid evaluation of the function f(grid)
       * @param x point to be interpolated
       * @param nearest_grid_index_to_x the nearest point to x is always
       * truncated with floor
       */
      inline double interpolate(const Vector_Ref & grid,
                                const Vector_Ref & evaluated_grid, double x,
                                int nearest_grid_index_to_x) {
        throw std::runtime_error(
            "Function `interpolate` for functions of the form f:[x1,x2]->ℝ has "
            "not been implemented in InterpolationMethod.");
        return 0;
      }

      /**
       * intp(x) for functions of the form f:[x1,x2]->ℝ^{n,m}.
       *
       * @param grid the grid for interpolation
       * @param evaluated_grid evaluation of the function f(grid) of shape
       *                       (grid_size, n*m)
       * @param x point to be interpolated
       * @param nearest_grid_index_to_x the nearest point to x is always
       * truncated with floor
       */
      inline double interpolate(const Vector_Ref & grid,
                                const Matrix_Ref & evaluated_grid, double x,
                                int nearest_grid_index_to_x) {
        throw std::runtime_error(
            "Function `interpolate` for functions of the form f:[x1,x2]->ℝ^{n,m} has "
            "not been implemented in InterpolationMethod.");
        return 0;
      }
    };

    /**
     * An implementation of the cubic spline method, the implementation is
     * adapted from "Numerical Recipes" [1] for functions of the form
     * f:[x1,x2]->ℝ and optimized for uniform grids, works only for uniform grids.
     *
     * [1] Press, William H., et al. Numerical recipes 3rd edition: The art of
     * scientific computing. Cambridge university press, 2007.
     */
    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSplineScalarUniform> {
     public:
      constexpr static InterpolationMethod_t Method{
          InterpolationMethod_t::CubicSplineScalarUniform};

      InterpolationMethod<InterpolationMethod_t::CubicSplineScalarUniform>() {}

      void initialize(const Vector_Ref & grid,
                      const Vector_Ref & evaluated_grid) {
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h * this->h / 6.0;
        this->compute_second_derivative(evaluated_grid);
      }

      void initialize(const Vector_Ref & grid,
                      const Vector_Ref & evaluated_grid, double yd1,
                      double ydn) {
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h * this->h / 6.0;
        this->compute_second_derivative(evaluated_grid, yd1, ydn);
      }

      inline double interpolate(const Vector_Ref & grid,
                                const Vector_Ref & evaluated_grid, double x,
                                int nearest_grid_index_to_x) {
        return this->raw_interpolate(grid, evaluated_grid,
                                     nearest_grid_index_to_x, x);
      }
      inline double interpolate_derivative(const Vector_Ref & grid,
                                           const Vector_Ref & evaluated_grid,
                                           double x,
                                           int nearest_grid_index_to_x) {
        return this->raw_interpolate_derivative(grid, evaluated_grid,
                                                nearest_grid_index_to_x, x);
      }

     private:
      /**
       * The two functions below implement the tridiagonal matrix algorithm
       * to solve the linear system in
       * http://mathworld.wolfram.com/CubicSpline.html
       * under different boundary conditions. We describe with s the spline
       * function and with y the targeted function for interpolation.
       *
       * Reference:
       * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

       * clamped boundary conditions when first derivative for boundaries is
       * known s'(x_0) = y'(x_0), s'(x_n) = y'(x_n)
       */
      inline void compute_second_derivative(const Vector_Ref & yv, double yd1,
                                            double ydn) {
        // Bad xa input to routine splint
        assert(this->h > dbl_ftol);
        int n{static_cast<int>(yv.size())};
        Vector_t y2 = Vector_t::Zero(n);
        Vector_t u = Vector_t::Zero(n);
        double p;

        // forward sweeping
        y2(0) = -0.5;
        u(0) = (3.0 / this->h) * ((yv(1) - yv(0)) / this->h - yd1);
        for (int i{1}; i < n - 1; i++) {
          p = 0.5 * y2(i - 1) + 2.0;
          y2(i) = -0.5 / p;
          u(i) = (yv(i + 1) - 2 * yv(i) + yv(i - 1)) / this->h;
          u(i) = (3.0 * u(i) / this->h - 0.5 * u(i - 1)) / p;
        }

        // back substitution
        p = 0.5;  // qn=0.5 but we just reuse p because it is not needed anyway
        u(n - 1) = (3.0 / this->h) * (ydn - (yv(n - 1) - yv(n - 2)) / this->h);
        y2(n - 1) = (u(n - 1) - p * u(n - 2)) / (p * y2(n - 2) + 1.0);
        for (int k{n - 2}; k >= 0; k--) {
          y2(k) = y2(k) * y2(k + 1) + u(k);
        }
        this->second_derivative_h_sq_6 = y2 * this->h_sq_6;
      }

      // natural/simple boundary conditions s''(x_0) = s''(x_n) = 0
      inline void compute_second_derivative(const Vector_Ref & yv) {
        // Bad xa input to routine splint
        assert(this->h > dbl_ftol);
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

      inline double raw_interpolate(const Vector_Ref & xx,
                                    const Vector_Ref & yy, int j1,
                                    double x) {
        assert(this->h > dbl_ftol);
        // j1
        const int klo{j1}, khi{j1 + 1};
        // percentage of grid cell start to x
        const double a{(xx(khi) - x) / this->h};
        // percentage of x to grid cell end, because a+b = 1 we can simplify
        // the computation
        double b{1 - a};
        return a * (yy(klo) +
                    (a * a - 1) * this->second_derivative_h_sq_6(klo)) +
               b * (yy(khi) +
                    (b * b - 1) * this->second_derivative_h_sq_6(khi));
      }

      inline double raw_interpolate_derivative(const Vector_Ref & xx,
                                               const Vector_Ref & yy,
                                               int j1,
                                               double x) {
        assert(this->h > dbl_ftol);
        const int klo{j1}, khi{j1 + 1};
        const double a{(xx(khi) - x) / this->h};
        // Because we can assume a+b=1, we simplify the calculation of b
        const double b{1 - a};
        return (yy(khi) - yy(klo) -
                (3 * a * a - 1) * this->second_derivative_h_sq_6(klo) +
                (3 * b * b - 1) * this->second_derivative_h_sq_6(khi)) /
               this->h;
      }

      // grid step size
      double h{0};
      // h*h/6.0
      double h_sq_6{0};
      /**
       * This term is the solution vector of the linear system of the
       * tridiagonal toeplitz matrix derived by the conditions of cubic spline.
       * The linear system can be found in
       * http://mathworld.wolfram.com/CubicSpline.html These terms are the
       * remaining constants of the second derivative of the cubic spline. In
       * addition we already multiply them with precomputed coefficients derived
       * by assuming a uniform grid y2 *h*h/6
       */
      Vector_t second_derivative_h_sq_6{};
    };

    /**
     * An implementation of the cubic spline method, the implementation is
     * adapted from "Numerical Recipes" [1] for functions of the form
     * y:[x1,x2]->ℝ^n and optimized for uniform grids.
     *
     * [1] Press, William H., et al. Numerical recipes 3rd edition: The art of
     * scientific computing. Cambridge university press, 2007.
     */
    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSplineVectorUniform> {
     public:
      constexpr static InterpolationMethod_t Method{
          InterpolationMethod_t::CubicSplineVectorUniform};

      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorUniform>() {}

      void initialize(const Vector_Ref & grid,
                      const Matrix_Ref & evaluated_grid) {
        this->h = grid(1) - grid(0);
        this->h_sq_6 = this->h * this->h / 6.0;
        this->compute_second_derivative(evaluated_grid);
      }

      inline Vector_t interpolate(const Vector_Ref & grid,
                                  const Matrix_Ref & evaluated_grid, double x,
                                  int nearest_grid_index_to_x) {
        return this->raw_interpolate(grid, evaluated_grid,
                                     nearest_grid_index_to_x, x);
      }
      inline Vector_t interpolate_derivative(const Vector_Ref & grid,
                                             const Matrix_Ref & evaluated_grid,
                                             double x,
                                             int nearest_grid_index_to_x) {
        return this->raw_interpolate_derivative(grid, evaluated_grid,
                                                nearest_grid_index_to_x, x);
      }

     private:
      inline void compute_second_derivative(const Matrix_Ref & yv) {
        int n{static_cast<int>(yv.rows())};
        Matrix_t y2 = Matrix_t::Zero(n, yv.cols());
        Matrix_t u = Matrix_t::Zero(n, yv.cols());
        const double sig = 0;  // wtf this is more accurate than the correct
                               // value 0.5 TODO(alex)
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
          y2.row(k) =
              y2.row(k).array() * y2.row(k + 1).array() + u.row(k).array();
        }
        y2.row(0) = Vector_t::Zero(y2.cols());
        this->second_derivative_h_sq_6 = y2 * this->h_sq_6;
      }

      inline Vector_t raw_interpolate(const Vector_Ref & xx,
                                      const Matrix_Ref & yy, int j1,
                                      double x) {
        int klo{j1}, khi{j1 + 1};
        // Bad xa input to routine splint
        assert(h != 0.0);
        double a{(xx(khi) - x) / this->h};
        // Because we can assume a+b=1, we simplify the calculation of b
        double b{1 - a};
        return (a * (yy.row(klo).array() +
                     (a * a - 1) *
                         this->second_derivative_h_sq_6.row(klo).array()) +
                b * (yy.row(khi).array() +
                     (b * b - 1) *
                         this->second_derivative_h_sq_6.row(khi).array()))
            .matrix();
      }

      inline Vector_t raw_interpolate_derivative(const Vector_Ref & xx,
                                                 const Matrix_Ref & yy,
                                                 int j1,
                                                 double x) {
        int klo{j1}, khi{j1 + 1};
        // Bad xa input to routine splint
        assert(h != 0.0);
        double a{(xx(khi) - x) / this->h};
        // Because we can assume a+b=1, we simplify the calculation of b
        double b{1 - a};
        return ((yy.row(khi).array() - yy.row(klo).array() -
                 (3 * a * a - 1) *
                     this->second_derivative_h_sq_6.row(klo).array() +
                 (3 * b * b - 1) *
                     this->second_derivative_h_sq_6.row(khi).array()) /
                this->h)
            .matrix();
      }

      // The gap between two grid points
      double h{0};
      // h*h/6.0
      double h_sq_6{0};
      // second_derivative*h*h/6
      Matrix_t second_derivative_h_sq_6{};
    };

    enum class ErrorMetric_t { Absolute, Relative };
    enum class ErrorNorm_t { Mean, Max };

    template <ErrorMetric_t Metric>
    struct ErrorMetric {
      /**
       * Computes the error of a function of the form f:[x1,x2]->ℝ of each point
       * on the grid.
       *
       * @param values intp(grid) by interpolator
       * @param references f(grid) by the given function
       */
      static Vector_t compute_entrywise_error(const Vector_Ref & values,
                                              const Vector_Ref & references) {
         throw std::runtime_error(
             "Function `compute_entrywise_error` for a function of the form "
             "f:[x1,x2]->ℝ has not been "
             "implemented in ErrorMetric.");
         return 0;
      }

      /**
       * Computes the error of a function of the form f:[x1,x2]->ℝ^{n,m} of each
       * point on the grid.
       *
       * @param values intp(grid) by interpolator
       * @param references f(grid) by the given function
       */
      static Matrix_t compute_entrywise_error(const Matrix_Ref & values,
                                              const Matrix_Ref & references) {
         throw std::runtime_error(
             "Function `compute_entrywise_error` for a function of the form "
             "f:[x1,x2]->ℝ^{n,m} has not been "
             "implemented in ErrorMetric.");
         return 0;
      }
    };

    template <ErrorMetric_t Metric, ErrorNorm_t Norm>
    struct ErrorMethod {
      /**
       * Computes the error of a function of the form f:[x1,x2]->ℝ of the whole grid.
       *
       * @param values intp(grid) by interpolator
       * @param references f(grid) by the given function
       */
      static double compute_global_error(const Vector_Ref & values, const
                                         Vector_Ref & references) {
        throw std::runtime_error(
            "Function `compute_global_error` for a function of the form "
            "f:[x1,x2]->ℝ has not been "
            "implemented in ErrorMethod.");
        return 0;
      }

      /**
       * Computes the error for a grid with a function of the form
       * f:[x1,x2]->ℝ^{n,m} reduced to a scalar value.
       *
       * @param values f(grid) by interpolator with shape (grid_size, n*m)
       * @param references f(grid) by the given function with shape
       *                   (grid_size, n*m)
       *
       */
      static double compute_global_error(const Matrix_Ref & values, const
      Matrix_Ref & references) {
        throw std::runtime_error(
            "Function `compute_global_error` for a function of the form "
            "f:[x1,x2]->ℝ{n,m} has not been "
            "implemented in ErrorMethod.");
        return 0;
      }
    };

    /**
     * Computes the relative error entrywise.
     *
     * Works well for interpolators with functions evaluating large absolute
     * values, outside the range [-1,1], so for functions of the form
     * f:[x1,x2]->[a,b] with a << -1 and b >> 1.
     */
    template <>
    struct ErrorMetric<ErrorMetric_t::Relative> {
      static Vector_t compute_entrywise_error(const Vector_Ref & values,
                                              const Vector_Ref & references) {
        return 2 * ((values - references).array().abs() /
                    (std::numeric_limits<double>::min() + values.array().abs() +
                     references.array().abs()));
      }

      static Matrix_t compute_entrywise_error(const Matrix_Ref & values,
                                              const Matrix_Ref & references) {
        return 2 * ((values - references).array().abs() /
                    (std::numeric_limits<double>::min() + values.array().abs() +
                     references.array().abs()));
      }
    };

    /**
     * Computes the absolute error entrywise.
     *
     * Works well for interpolators with functions evaluating small absolute
     * values, inside the range [-1,1], so for functions of the form
     * f:[x1,x2]->[-1,1].
     */
    template <>
    struct ErrorMetric<ErrorMetric_t::Absolute> {
      static Vector_t compute_entrywise_error(const Vector_Ref & values,
                                              const Vector_Ref & references) {
        return (values - references).array().abs();
      }

      static Matrix_t compute_entrywise_error(const Matrix_Ref & values,
                                              const Matrix_Ref & references) {
        return (values - references).array().abs();
      }
    };

    /**
     * Computes the mean error. For functions of the form f:[x1,x2]->ℝ^{n,m} the
     * maximum of the n*m mean errors is taken.
     */
    template <ErrorMetric_t Metric>
    struct ErrorMethod<Metric, ErrorNorm_t::Mean> : public ErrorMetric<Metric> {
      using Parent = ErrorMetric<Metric>;
      static double compute_global_error(const Vector_Ref & values,
                                         const Vector_Ref & references) {
        return Parent::compute_entrywise_error(values, references).mean();
      }
      static double compute_global_error(const Matrix_Ref & values,
                                         const Matrix_Ref & references) {
        return Parent::compute_entrywise_error(values, references)
            .colwise()
            .mean()
            .maxCoeff();
      }
    };

    template <ErrorMetric_t Metric>
    struct ErrorMethod<Metric, ErrorNorm_t::Max> : public ErrorMetric<Metric> {
      using Parent = ErrorMetric<Metric>;
      static double compute_global_error(const Vector_Ref & values,
                                         const Vector_Ref & references) {
        return Parent::compute_entrywise_error(values, references).maxCoeff();
      }

      static double compute_global_error(const Matrix_Ref & values,
                                         const Matrix_Ref & references) {
        return Parent::compute_entrywise_error(values, references).maxCoeff();
      }
    };

    enum class SearchMethod_t { Hunt, Locate, Uniform };

    template <SearchMethod_t Method>
    struct SearchMethod {
      /**
       * The main purpose of the Search method is to obtain from a point x in
       * the range [x1,x2] the nearest grid point. With nearest the grid point
       * just before point x is meant. x1     x        x2
       *                            |
       * grid                [   |   |   |   ]
       * nearest grid point      |
       *
       * @param grid used for searching the index
       * void initialize(const Vector_Ref & grid)
       *
       * @param grid used for searching the index
       * @param x as mentioned above
       * @return nearest grid point
       * int search(double x, const Vector_Ref & grid)
       */
    };

    template <>
    class SearchMethod<SearchMethod_t::Uniform> {
     public:
      constexpr static SearchMethod_t Method{SearchMethod_t::Uniform};

      void initialize(const Vector_Ref & grid) {
        // nb_grid_points/unit
        this->nb_grid_points_per_unit =
            grid.size() / (grid(grid.size() - 1) - grid(0));
        this->x1 = grid(0);
        this->grid_size = grid.size();
        this->search_size = this->grid_size - this->nb_support_points;
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref &) {
        return std::max(
            0, std::min(this->search_size,
                        static_cast<int>((x - this->x1) *
                                         this->nb_grid_points_per_unit) -
                            1));
      }

     private:
      // the number of support methods the interpolation method uses
      double nb_grid_points_per_unit{0};
      double x1{0};
      // asserts that we only use it with CubicSpline
      size_t nb_support_points{2};
      size_t grid_size{0};
      int search_size{0};
    };

    // Improvements:
    // - make nb_support_points const
    template <>
    struct SearchMethod<SearchMethod_t::Locate> {
      constexpr static SearchMethod_t Method{SearchMethod_t::Locate};

      SearchMethod<SearchMethod_t::Locate>() : nb_support_points{2} {}

      void initialize(const Vector_Ref &) {}

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        return this->locate(x, grid);
      }

      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};

        int ju, jm, jl;
        // size error
        assert(not(n < 2 || mm < 2 || mm > n));
        bool ascnd = (xx[n - 1] >= xx[0]);
        jl = 0;
        ju = n - 1;
        while (ju - jl > 1) {
          jm = (ju + jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl = jm;
          else
            ju = jm;
        }
        return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
      }

      // the number of support methods the interpolation method uses
      size_t nb_support_points;
    };

    // Improvements/Problems of hunt:
    // - Hunt gives indices outside of the range of the grid when compiled in
    // Debug or RelWithDebugInfo
    // - make nb_support_points const
    template <>
    struct SearchMethod<SearchMethod_t::Hunt> {
      constexpr static SearchMethod_t Method{SearchMethod_t::Hunt};
      SearchMethod<SearchMethod_t::Hunt>()
          : correlated{false}, nb_support_points{2},
            last_accessed_index{0}, dj{0} {}

      void initialize(const Vector_Ref & grid) {
        this->dj = std::min(
            1, static_cast<int>(std::round(std::sqrt(std::sqrt(grid.size())))));
      }

      // If the requests to locate seem correlated, then the heuristic is used
      int search(double x, const Vector_Ref & grid) {
        return this->correlated ? this->hunt(x, grid) : this->locate(x, grid);
      }

      int locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        int ju, jm, jl;
        assert(not(n < 2 || mm < 2 || mm > n));
        bool ascnd = (xx[n - 1] >= xx[0]);
        jl = 0;
        ju = n - 1;
        while (ju - jl > 1) {
          jm = (ju + jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl = jm;
          else
            ju = jm;
        }
        this->correlated = abs(jl - jsav) > this->dj ? 0 : 1;
        jsav = jl;

        this->last_accessed_index = jsav;
        return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
      }

      int hunt(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        int jl = jsav, jm, ju, inc = 1;
        // hunt size error
        assert(not(n < 2 || mm < 2 || mm > n));

        bool ascnd = (xx[n - 1] >= xx[0]);
        if (jl < 0 || jl > n - 1) {
          jl = 0;
          ju = n - 1;
        } else {
          if ((x >= xx[jl]) == ascnd) {
            for (;;) {
              ju = jl + inc;
              if (ju >= n - 1) {
                ju = n - 1;
                break;
              } else if ((x < xx[ju]) == ascnd) {
                break;
              } else {
                jl = ju;
                inc += inc;
              }
            }
          } else {
            ju = jl;
            for (;;) {
              jl = jl - inc;
              if (jl <= 0) {
                jl = 0;
                break;
              } else if ((x >= xx[jl]) == ascnd) {
                break;
              } else {
                ju = jl;
                inc += inc;
              }
            }
          }
        }
        while (ju - jl > 1) {
          jm = (ju + jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl = jm;
          else
            ju = jm;
        }
        this->correlated = abs(jl - jsav) > this->dj ? 0 : 1;
        jsav = jl;
        this->last_accessed_index = jsav;
        return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
      }

      bool correlated;
      // the number of support methods the interpolation method uses
      size_t nb_support_points;
      size_t last_accessed_index;
      // parameter used to determine if search requests are correlated
      int dj;
    };

    /**
     * Base class of an interpolator to allow storage of variable templates of
     * an interpolator as member variable, because templation of member
     * variables is not allowed. The interpolator upcasted for storage and
     * downcasted when used.
     */
    class InterpolatorBase {
     public:
      //! Constructor
      InterpolatorBase() = default;

      //! Destructor
      virtual ~InterpolatorBase() = default;

      //! Copy constructor
      InterpolatorBase(const InterpolatorBase & other) = delete;

      //! Move constructor
      InterpolatorBase(InterpolatorBase && other) = default;

      //! Copy assignment operator
      InterpolatorBase & operator=(const InterpolatorBase & other) = delete;

      //! Move assignment operator
      InterpolatorBase & operator=(InterpolatorBase && other) = default;
    };

    /**
     * Templated interpolator class, used as basic class for different kind of
     * interpolators (currently scalar, vector) to reduce code repetition and
     * allow additional optimization for specific types of interpolators.
     */
    template <class InterpolationMethod, class GridRational, class SearchMethod,
              class ErrorMethod>
    class Interpolator : public InterpolatorBase {
     public:
      /**
       * Different constructors for the interpolator, to be able to run the
       * interpolator with only the required parameters, but also to fine tune
       * different functionalities for an expert use. The detailed explanation
       * of each parameter can be seen in the protected constructor. Sets the
       * interpolator interpolating a function of the form y:[x1,x2]->Ω up to an
       * error bound or a max number of grid points in the range [x1,x2], where
       * Ω depends on the child inheriting from this class. Currently, there are
       * implementations for Ω equal ℝ or ℝ^n.
       *
       * @param function the function to be interpolated referred as y here
       * @param x1 begin range
       * @param x2 end range
       * @param error_bound the minimal error bound fulfilled
       * @param max_grid_points maximal number of grid points for the
       * interpolating grid until the refinement of the grid stops.
       * @param start_degree_of_fineness starting fineness of the grid when starting
       * the algorithm, the fineness parameter describes the fineness of the
       * grid for grid rationals with non adaptive refinement methods. Higher
       * value results in finer grids. It depends on the specifics of the
       * GridRational how the fineness changes the number of grid points
       * exactly. For adaptive refinement methods the fineness is ignored. The
       * fineness parameter should be at least 1. A good fineness value can
       * reduce the number of steps in the initialization.
       * @param empty bool parameter to prevent overloading with the public
       * constructor
       *
       * @error range is illdefined
       * @error error bound is illdefined
       * @error fineness is illdefined.
       */
      Interpolator(double x1, double x2, double error_bound,
                   int max_grid_points, int start_degree_of_fineness)
          : Interpolator(x1, x2, error_bound, max_grid_points,
                         start_degree_of_fineness, true) {}
      Interpolator(double x1, double x2, double error_bound)
          : Interpolator(x1, x2, error_bound, 10000000, 5, true) {}

      /**
       * Constructor to initialize interpolator from grid.
       * @param grid in the range [x1,x2]
       */
      explicit Interpolator(Vector_t grid) : Interpolator(grid, true) {}

      int get_degree_of_fineness() { return this->degree_of_fineness; }
      int get_grid_size() { return this->grid.size(); }
      Vector_Ref get_grid_ref() { return Vector_Ref(this->grid); }

      // Should be implemented, but is not made virtual for performance
      // virtual Ω interpolate(double x);
      // virtual Ω interpolate_derivative(double x);
      // virtual Vector_Ref get_grid_ref() = 0;
      // virtual Ω  get_evaluated_grid_ref() = 0;
      // virtual int get_grid_size() = 0;

     protected:
      Interpolator(double x1, double x2, double error_bound,
                   int max_grid_points, int start_degree_of_fineness, bool)
          : x1{x1}, x2{x2}, error_bound{error_bound},
            max_grid_points{max_grid_points},
            degree_of_fineness{start_degree_of_fineness},
            intp_method{InterpolationMethod()}, grid_rational{GridRational()},
            search_method{SearchMethod()} {
        if (x2 < x1) {
          throw std::logic_error("x2 must be greater x1");
        }
        if (error_bound <= 0) {
          throw std::logic_error("Error bound must be > 0");
        }
        if (start_degree_of_fineness < 1) {
          throw std::logic_error("Starting grid fineness must be at least 1.");
        }
        if (max_grid_points < 2) {
          throw std::logic_error(
              "Maximal number of grid points must be at least 2.");
        }
      };

      Interpolator(Vector_t grid, bool)
          : grid{grid}, intp_method{InterpolationMethod()},
            grid_rational{GridRational()}, search_method{SearchMethod()} {}

      /**
       * The general procedure of the interpolatorar is to initialize the
       * ressources for the interpolation method on a grid and test if the error
       * is within a bound. If this is not the case the interpolator will refine
       * the grid and repeat its procedure.
       */
      void initialize_iteratively() {
        this->compute_grid_error();
        while (this->grid_error > this->error_bound &&
               this->grid.size() < this->max_grid_points) {
          this->degree_of_fineness++;
          this->compute_grid_error();
        }
      }

      virtual void compute_grid_error() = 0;

      // The boundary points of the range of interpolation
      const double x1{0.};
      const double x2{0.};
      // grid in the range [x1,x2]
      Vector_t grid{};
      const double error_bound{1e-5};
      double grid_error{0.};

      const int max_grid_points{10000000};  // 1e7
      int degree_of_fineness{5};

      InterpolationMethod intp_method{};
      GridRational grid_rational{};
      SearchMethod search_method{};
    };

    /**
     * Specialization of the interpolator for uniform grids with cubic spline.
     * The spline methods can be also adapted for general grids, however for
     * uniform grids the computation can be much more optimized.
     */
    template <RefinementMethod_t RefinementMethod,
              class ErrorMethod_ =
                  ErrorMethod<ErrorMetric_t::Absolute, ErrorNorm_t::Mean>>
    class InterpolatorScalarUniformCubicSpline
        : public Interpolator<
              InterpolationMethod<
                  InterpolationMethod_t::CubicSplineScalarUniform>,
              GridRational<RefinementMethod>,
              SearchMethod<SearchMethod_t::Uniform>, ErrorMethod_> {
     public:
      using Parent = Interpolator<
          InterpolationMethod<InterpolationMethod_t::CubicSplineScalarUniform>,
          GridRational<RefinementMethod>,
          SearchMethod<SearchMethod_t::Uniform>, ErrorMethod_>;
      using This =
          InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod_>;
      using ThisErrorMethod = ErrorMethod_;
      using ThisGridRational =
          GridRational<RefinementMethod>;
      using ThisSearchMethod = SearchMethod<SearchMethod_t::Uniform>;

      /**
       * Constructor using a the function of the form y:[x1,x2]->ℝ to initialize
       * a grid by iteratively refining it until an error bound (logical) or a
       * max number of points for the grid is reached.
       *
       * @param function the function to be interpolated referred as y here
       * @param x1 begin range
       * @param x2 end range
       * @param error_bound the minimal error bound fulfilled
       * @param max_grid_points maximal number of grid points for the
       * interpolating grid until the refinement of the grid stops.
       * @param start_degree_of_fineness starting fineness of the grid when starting
       * the algorithm, the fineness parameter describes the fineness of the
       * grid for grid rationals with non adaptive refinement methods. Higher
       * value results in finer grids. It depends on the specifics of the
       * GridRational how the fineness changes the number of grid points
       * exactly. For adaptive refinement methods the fineness is ignored. The
       * fineness parameter should be at least 1. A good fineness value can
       * reduce the number of steps in the initialization.
       * @param clamped_boundary_conditions Parameter referring to the boundary
       * condition of cubic spline. By default the cubic spline uses the natural
       * boundary conditions, but if y'(x1) and y'(x2) are known, the clamped
       * boundary conditions can be turned on, which can increase the accuracy
       * of the interpolation.
       * @param yd1 referring to y'(x1)
       * @param ydn referring to y'(x2)
       */
      InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          std::function<double(double)> function, double x1, double x2,
          double error_bound)
          : Parent(x1, x2, error_bound), function{function},
            clamped_boundary_conditions{false}, yd1{0}, ydn{0} {
        this->initialize_iteratively();
      }

      InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          std::function<double(double)> function, double x1, double x2,
          double error_bound, int max_grid_points, int start_degree_of_fineness,
          bool clamped_boundary_conditions, double yd1, double ydn)
          : Parent{x1, x2, error_bound, max_grid_points, start_degree_of_fineness},
            function{function},
            clamped_boundary_conditions{clamped_boundary_conditions}, yd1{yd1},
            ydn{ydn} {
        this->initialize_iteratively();
      }

      /**
       * Constructor using a grid to initialize the interpolator in one step
       * interpolating a function of the form y:[x1,x2]->ℝ.
       *
       * @param grid a uniform grid in the range [x1,x2] with x1 and x2 as start
       * @param evaluated_grid the grid evaluated by the function referred as
       * y(grid)
       * @param clamped_boundary_conditions By default the cubic spline uses the
       * natural boundary conditions, but if y'(x1) and y'(x2) are known, the
       * clamped boundary conditions can be turned on, which can increase the
       * accuracy of the interpolation.
       */
      InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          Vector_t grid, Vector_t evaluated_grid)
          : Parent(grid), evaluated_grid{evaluated_grid},
            clamped_boundary_conditions{false}, yd1{0}, ydn{0} {
        if (grid.size() != evaluated_grid.size()) {
          throw std::logic_error(
              "The grid and evaluated grid must match in size");
        }
        if (not(this->grid_rational.is_grid_uniform(Vector_Ref(this->grid)))) {
          throw std::logic_error("The grid has to be unform.");
        }
        this->initialize_from_computed_grid();
      }

      InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          Vector_t grid, Vector_t evaluated_grid,
          bool clamped_boundary_conditions, double yd1, double ydn)
          : Parent(grid), evaluated_grid{evaluated_grid},
            clamped_boundary_conditions{clamped_boundary_conditions}, yd1{yd1},
            ydn{ydn} {
        if (grid.size() != evaluated_grid.size()) {
          throw std::logic_error(
              "The grid and evaluated grid must match in size");
        }
        if (not(this->grid_rational.is_grid_uniform(Vector_Ref(this->grid)))) {
          throw std::logic_error("The grid has to be unform.");
        }
        this->initialize_from_computed_grid();
      }

      double interpolate(double x) {
        assert(x >= this->x1 && x <= this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate(this->grid, this->evaluated_grid,
                                             x, nearest_grid_index_to_x);
      }

      Vector_t interpolate(const Vector_Ref & points) {
        Vector_t interpolated_points = Vector_t::Zero(points.size());
        for (int i{0}; i < points.size(); i++) {
          interpolated_points(i) = this->interpolate(points(i));
        }
        return interpolated_points;
      }

      // OPT(alex) search method vectorize
      double interpolate_derivative(double x) {
        assert(x >= this->x1 && x <= this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate_derivative(
            this->grid, this->evaluated_grid, x, nearest_grid_index_to_x);
      }

      Vector_Ref get_evaluated_grid_ref() {
        return Vector_Ref(this->evaluated_grid);
      }

      /**
       * Function returns a json with information about the interpolator to
       * evaluate the interpolator. Used for benchmarks.
       *
       * @return json with information about the interpolator
       */
      json compute_interpolator_information() {
        // compute test errors
        Vector_t test_grid{this->grid_rational.compute_test_grid(
            this->x1, this->x2, this->degree_of_fineness)};
        Vector_t test_grid_interpolated{
            this->interpolate(Vector_Ref(test_grid))};
        Vector_t test_grid_evaluated{this->eval(Vector_Ref(test_grid))};
        Vector_t error_grid{ThisErrorMethod::compute_entrywise_error(
            Vector_Ref(test_grid_interpolated),
            Vector_Ref(test_grid_evaluated))};
        double mean_grid_error{error_grid.mean()};
        double max_grid_error{error_grid.maxCoeff()};

        return {{"mean_grid_error", mean_grid_error},
                {"max_grid_error", max_grid_error},
                {"grid_size", this->grid.size()}};
      }

     protected:
      void initialize_from_computed_grid() {
        if (this->clamped_boundary_conditions) {
          this->intp_method.initialize(Vector_Ref(this->grid),
                                       Vector_Ref(this->evaluated_grid),
                                       this->yd1, this->ydn);
        } else {
          this->intp_method.initialize(Vector_Ref(this->grid),
                                       Vector_Ref(this->evaluated_grid));
        }
        this->search_method.initialize(Vector_Ref(this->grid));
      }

      /**
       * Estimates the interpolator error by computing a refined test grid and
       * compares the result of the function with the interpolator
       */
      void compute_grid_error() {
        // compute the grid
        this->grid = this->grid_rational.compute_grid(this->x1, this->x2,
                                                      this->degree_of_fineness);
        // OPT(alex) use the test grid to save computation time in evaluating
        // the grid
        this->evaluated_grid = this->eval(Vector_Ref(this->grid));
        this->initialize_from_computed_grid();

        // compute test grid
        Vector_t test_grid{this->grid_rational.compute_test_grid(
            this->x1, this->x2, this->degree_of_fineness)};
        Vector_t test_grid_interpolated{
            this->interpolate(Vector_Ref(test_grid))};
        Vector_t test_grid_evaluated{this->eval(Vector_Ref(test_grid))};

        // error method ε:ℝ^{grid_points}->ℝ
        this->grid_error = ThisErrorMethod::compute_global_error(
            Vector_Ref(test_grid_interpolated),
            Vector_Ref(test_grid_evaluated));

      }

      double eval(double x) { return this->function(x); }

      // We use evaluate when the function is used and interpolate when the
      // interpolation method is used
      Vector_t eval(const Vector_Ref & grid) {
        Vector_t evaluated_grid = Vector_t::Zero(grid.size());
        for (int i{0}; i < evaluated_grid.size(); i++) {
          evaluated_grid(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }

      /**
       * This computes parameters which are needed in the benchmarks for
       * evaluation.
       */
      void compute_paratemeters_for_evaluation() {
        Vector_t test_grid{this->grid_rational.compute_test_grid(
            this->x1, this->x2, this->degree_of_fineness)};
        Vector_t test_grid_interpolated{this->interpolate(test_grid)};
        Vector_t test_grid_evaluated{this->eval(test_grid)};
        Vector_t error_grid{ThisErrorMethod::compute_entrywise_error(
            Vector_Ref(test_grid_interpolated),
            Vector_Ref(test_grid_evaluated))};
        this->max_grid_error = error_grid.maxCoeff();
        this->mean_grid_error = error_grid.mean();
      }
      // f:[x1,x2]->ℝ
      std::function<double(double)> function{};
      // evaluated by the function to interpolate f(grid)
      Vector_t evaluated_grid{};
      double mean_grid_error{0.};
      double max_grid_error{0.};
      /**
       * If the first derivative of the function at the boundary points is
       * known, one can initialize the interpolator with these values to obtain
       * a higher accuracy with the interpolation method.
       */
      const bool clamped_boundary_conditions{false};
      // f'(x1)
      const double yd1{0.};
      // f'(x2)
      const double ydn{0.};
    };

    template <RefinementMethod_t RefinementMethod,
              class ErrorMethod_ =
                  ErrorMethod<ErrorMetric_t::Absolute, ErrorNorm_t::Mean>>
    class InterpolatorMatrixUniformCubicSpline
        : public Interpolator<
              InterpolationMethod<
                  InterpolationMethod_t::CubicSplineVectorUniform>,
              GridRational<RefinementMethod>,
              SearchMethod<SearchMethod_t::Uniform>, ErrorMethod_> {
     public:
      using Parent = Interpolator<
          InterpolationMethod<InterpolationMethod_t::CubicSplineVectorUniform>,
          GridRational<RefinementMethod>,
          SearchMethod<SearchMethod_t::Uniform>, ErrorMethod_>;
      using ThisErrorMethod = ErrorMethod_;
      using ThisGridRational =
          GridRational<RefinementMethod>;

      /**
       * Constructor using a the function of the form y:[x1,x2]->ℝ^{n,m} to
       * initialize a grid by iteratively refining it until an error bound
       * (logical) or a max number of points for the grid is reached.
       *
       * @param function the function to be interpolated referred as y here
       * @param x1 begin range
       * @param x2 end range
       * @param error_bound the minimal error bound fulfilled
       * @param max_grid_points maximal number of grid points for the
       * interpolating grid until the refinement of the grid stops.
       * @param start_degree_of_fineness starting fineness of the grid when starting
       * the algorithm, the fineness parameter describes the fineness of the
       * grid for grid rationals with non adaptive refinement methods. Higher
       * value results in finer grids. It depends on the specifics of the
       * GridRational how the fineness changes the number of grid points
       * exactly. For adaptive refinement methods the fineness is ignored. The
       * fineness parameter should be at least 1. A good fineness value can
       * reduce the number of steps in the initialization.
       * @param clamped_boundary_conditions Parameter referring to the boundary
       * condition of cubic spline. By default the cubic spline uses the natural
       * boundary conditions, but if y'(x1) and y'(x2) are known, the clamped
       * boundary conditions can be turned on, which can increase the accuracy
       * of the interpolation.
       * @param yd1 referring to y'(x1)
       * @param ydn referring to y'(x2)
       */
      InterpolatorMatrixUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          std::function<Matrix_t(double)> function, double x1,
          double x2, double error_bound, int cols,
          int rows)
          : Parent(x1, x2, error_bound), function{function}, cols{cols},
            rows{rows}, matrix_size{cols * rows},
            clamped_boundary_conditions{false}, yd1{0}, ydn{0} {
        this->initialize_iteratively();
      }

      InterpolatorMatrixUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          std::function<Matrix_t(double)> function, double x1, double x2,
          double error_bound, int cols, int rows, int max_grid_points,
          int start_degree_of_fineness, bool clamped_boundary_conditions, double yd1,
          double ydn)
          : Parent{x1, x2, error_bound, max_grid_points, start_degree_of_fineness},
            function{function}, cols{cols}, rows{rows}, matrix_size{cols *
                                                                    rows},
            clamped_boundary_conditions{clamped_boundary_conditions}, yd1{yd1},
            ydn{ydn} {
        if (clamped_boundary_conditions) {
          throw std::logic_error("Clamped boundary condition has yet not been "
                                 "implemented for CubicSplineVectorUniform.");
        }
        this->initialize_iteratively();
      }

      /**
       * Constructor using a grid to initialize the interpolator in one step
       * interpolating a function of the form y:[x1,x2]->ℝ^{cols, rows}
       *
       * @param grid a uniform grid in the range [x1,x2] with x1 and x2 as start
       * @param evaluated_grid the grid evaluated by the function referred as
       * y(grid) with shape (grid_size, cols*rows)
       * @param cols
       * @param rows
       * @param clamped_boundary_conditions By default the cubic spline uses the
       * natural boundary conditions, but if y'(x1) and y'(x2) are known, the
       * clamped boundary conditions can be turned on, which can increase the
       * accuracy of the interpolation.
       * @param yd1 referring to y'(x1)
       * @param ydn referring to y'(x2)
       *
       * @error grid is not uniform
       * @error grid and evaluated grid do not match by size
       * @error evaluated grid and cols and rows do not match
       */
      InterpolatorMatrixUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          Vector_t grid, Matrix_t evaluated_grid, int cols, int rows)
          : Parent(grid), evaluated_grid{evaluated_grid}, cols{cols},
            rows{rows}, matrix_size{cols * rows},
            clamped_boundary_conditions{false}, yd1{0}, ydn{0} {
        if (grid.size() != evaluated_grid.rows()) {
          throw std::logic_error(
              "The grid size and evaluated grid rows must match");
        }
        if (not(this->grid_rational.is_grid_uniform(Vector_Ref(this->grid)))) {
          throw std::logic_error("The grid has to be uniform.");
        }
        if (not(evaluated_grid.cols() == cols * rows)) {
          throw std::logic_error(
              "The evaluated grid number of cols must match cols*rows");
        }
        this->initialize_from_computed_grid();
      }

      InterpolatorMatrixUniformCubicSpline<RefinementMethod, ErrorMethod_>(
          Vector_t grid, Matrix_t evaluated_grid, int cols, int rows,
          bool clamped_boundary_conditions, double yd1, double ydn)
          : Parent(grid), evaluated_grid{evaluated_grid}, cols{cols},
            rows{rows}, matrix_size{cols * rows},
            clamped_boundary_conditions{clamped_boundary_conditions}, yd1{yd1},
            ydn{ydn} {
        if (clamped_boundary_conditions) {
          throw std::logic_error("Clamped boundary condition has yet not been "
                                 "implemented for CubicSplineVectorUniform.");
        }
        if (grid.size() != evaluated_grid.rows()) {
          throw std::logic_error(
              "The grid size and evaluated grid rows must match");
        }
        if (not(this->grid_rational.is_grid_uniform(Vector_Ref(this->grid)))) {
          throw std::logic_error("The grid has to be uniform.");
        }
        if (not(evaluated_grid.cols() == cols * rows)) {
          throw std::logic_error(
              "The evaluated grid number of cols must match cols*rows");
        }
        this->initialize_from_computed_grid();
      }

      // OPT(alex) container for Matrix_t, then reshape one time to prevent
      // check overflow
      /**
       * @error(debug) point is not in range [x1,x2]
       * @param x point to be interpolated within range [x1,x2]
       * @return interpolation of point x, estimating y(x)
       */
      Matrix_t interpolate(double x) {
        return Eigen::Map<Matrix_t>(this->raw_interpolate(x).data(), this->rows,
                                    this->cols);
      }

      /**
       * interpolates the derivative
       *
       * @error(debug) point is not in range [x1,x2]
       * @param x point to be interpolated within range [x1,x2]
       * @return interpolation of point x, estimating y'(x)
       */
      Matrix_t interpolate_derivative(double x) {
        return Eigen::Map<Matrix_t>(this->raw_interpolate_derivative(x).data(),
                                    this->rows, this->cols);
      }

      /**
       * This function exists for benchmark purposes, to estimate the memory
       * and function call overhead.
       *
       * @error(debug) point is not in range [x1,x2]
       * @param x point to be interpolated within range [x1,x2]
       * @return interpolation of point x, estimating y(x)
       */
      Matrix_t interpolate_dummy(const double x) const {
        assert(x >= this->x1 && x <= this->x2);
        return x * Matrix_t::Ones(this->rows, this->cols);
      }

      Matrix_Ref get_evaluated_grid_ref() {
        return Matrix_Ref(this->evaluated_grid);
      }
      int get_matrix_size() { return this->matrix_size; }

      /**
       * Function returns a json with information about the interpolator to
       * evaluate the interpolator. Used for benchmarks.
       *
       * @return json with information about the interpolator
       */
      json compute_interpolator_information() {
        // compute test errors
        Vector_t test_grid{this->grid_rational.compute_test_grid(
            this->x1, this->x2, this->degree_of_fineness)};
        Matrix_t test_grid_interpolated{
            this->raw_interpolate(Vector_Ref(test_grid))};
        Matrix_t test_grid_evaluated{this->eval(Vector_Ref(test_grid))};
        Matrix_t error_mat{ThisErrorMethod::compute_entrywise_error(
            Matrix_Ref(test_grid_interpolated),
            Matrix_Ref(test_grid_evaluated))};
        double mean_grid_error{error_mat.mean()};
        double max_grid_error{error_mat.maxCoeff()};

        return {{"mean_grid_error", mean_grid_error},
                {"max_grid_error", max_grid_error},
                {"grid_size", this->grid.size()}};
      }

     protected:
      // ensures that the grid for estimating the test error is fine enough.
      void compute_grid_error() {
        this->grid = this->grid_rational.compute_grid(this->x1, this->x2,
                                                      this->degree_of_fineness);
        this->evaluated_grid = this->eval(Vector_Ref(this->grid));
        this->initialize_from_computed_grid();

        // ensures that the grid for estimating the test error is fine enough.
        int test_degree_of_fineness{std::max(this->degree_of_fineness, 5)};
        Vector_t test_grid{this->grid_rational.compute_test_grid(
            this->x1, this->x2, test_degree_of_fineness)};
        // (grid_size, row*col)
        // OPT(alex) check how this is optimized
        Matrix_t test_grid_interpolated{
            this->raw_interpolate(Vector_Ref(test_grid))};
        Matrix_t test_grid_evaluated{this->eval(Vector_Ref(test_grid))};
        this->grid_error = ThisErrorMethod::compute_global_error(
            Matrix_Ref(test_grid_interpolated),
            Matrix_Ref(test_grid_evaluated));
      }

      void initialize_from_computed_grid() {
        if (this->clamped_boundary_conditions) {
          // OPT(alex) implement for clamped_boundary_conditions
          // this->intp_method.initialize(Vector_Ref(this->grid),
          // Matrix_Ref(this->evaluated_grid), this->yd1, this->ydn);
          this->intp_method.initialize(Vector_Ref(this->grid),
                                       Matrix_Ref(this->evaluated_grid));
        } else {
          this->intp_method.initialize(Vector_Ref(this->grid),
                                       Matrix_Ref(this->evaluated_grid));
        }
        this->search_method.initialize(Vector_Ref(this->grid));
      }

      Vector_t eval(double x) {
        return Eigen::Map<Vector_t>(this->function(x).data(),
                                    this->matrix_size);
      }

      // OPT(alex) container for Matrix_t, then reshape one time to prevent
      // check overflow
      Matrix_t eval(const Vector_Ref & grid) {
        Matrix_t evaluated_grid =
            Matrix_t::Zero(grid.size(), this->matrix_size);
        for (int i{0}; i < grid.size(); i++) {
          evaluated_grid.row(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }

      /**
       * Interpolates the point x
       * @return intp(x) in the shape (rows*cols)
       */
      inline Vector_t raw_interpolate(double x) {
        // x is outside of range
        assert(x >= this->x1 && x <= this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate(this->grid, this->evaluated_grid,
                                             x, nearest_grid_index_to_x);
      }

      /**
       * Interpolates the each x in points
       * @return intp(x) in the shape (points.size(), rows*cols)
       */
      inline Matrix_t raw_interpolate(const Vector_Ref & points) {
        Matrix_t interpolated_points =
            Matrix_t::Zero(points.size(), this->matrix_size);
        for (int i{0}; i < points.size(); i++) {
          interpolated_points.row(i) = this->raw_interpolate(points(i));
        }
        return interpolated_points;
      }

      Vector_t raw_interpolate_derivative(double x) {
        // x is outside of range
        assert(x >= this->x1 && x <= this->x2);
        int nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return this->intp_method.interpolate_derivative(
            this->grid, this->evaluated_grid, x, nearest_grid_index_to_x);
      }

      // TODO(alex) bring get information to Handler and Manager

      // y:[x1,x2]->ℝ^{rows,cols}
      std::function<Matrix_t(double)> function{};
      // y(grid) reshaped to (grid_size, n*m)
      Matrix_t evaluated_grid{};
      const int cols{0};
      const int rows{0};
      // cols*rows
      const int matrix_size{0};
      /**
       * If the first derivative of the function at the boundary points is
       * known, one can initialize the interpolator with these values to obtain
       * a higher accuracy with the interpolation method.
       */
      const bool clamped_boundary_conditions{false};
      // y'(x1)
      const double yd1{0.};
      // y'(x2)
      const double ydn{0.};
    };

  }  // namespace math
}  // namespace rascal
#endif  // SRC_MATH_INTERPOLATOR_HH_
