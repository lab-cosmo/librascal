/**
 * @file  performance/profiles/profile_scalar_cubic_spline.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief used to profile the scalar cubic spline method by storing the
 * grid in the initialization process in a file and loading it on start, to
 * remove the computation cost of the initialization process during the
 * profiling
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

#include "utils.hh"

#include "rascal/json_io.hh"
#include "rascal/math/hyp1f1.hh"
#include "rascal/math/interpolator.hh"

#include <chrono>
#include <iostream>

static constexpr int N_REPETITIONS = 200;
static unsigned int SEED = 1597463007;

using namespace rascal;  // NOLINT

using math::Hyp1f1;
using math::InterpolatorScalarUniformCubicSpline;
using math::RefinementMethod_t;
using math::Vector_t;

/**
 * If one wants to profile only the interpolate function, please execute this
 * file one time before using the profiler to create the grid with the
 * interpolator. The computed grid is then used in the second run. If one wants
 * to profile the interpolate function together with the grid initialization
 * process, be sure to delete the *.grid file before.
 */
int main() {
  // hyp1f1 parameters
  double n = 5;
  double l = 4;
  double a = 0.5 * (n + l + 3);
  double b = l + 1.5;
  auto hyp1f1 = Hyp1f1(a, b, 200, 1e-15);
  std::function<double(double)> func = [&hyp1f1](double x) {
    return hyp1f1.calc(x, 0.25, 0.5, true);
  };

  // interpolator parameters
  using IntpScalarUniformCubicSpline =
      InterpolatorScalarUniformCubicSpline<RefinementMethod_t::Exponential>;
  std::shared_ptr<IntpScalarUniformCubicSpline> intp;
  double x1{0};
  double x2{8};
  double error_bound{1e-5};

  // profile parameters
  size_t nb_points = 1e6;
  size_t nb_iterations = 1000000;
  const char * filename_grid{"profile_scalar_cubic_spline_grid.grid"};
  const char * filename_evaluated_grid{
      "profile_scalar_cubic_spline_grid.evaluated_grid"};

  // shuffle points to test interpolation method for uncorrelated requests
  srand(SEED);
  Vector_t points_tmp = Vector_t::LinSpaced(nb_points, x1, x2);
  Vector_t points = Vector_t::Zero(nb_points);
  for (size_t i{0}; i < nb_points; i++) {
    points(i) = points_tmp(rand_r(&SEED) % nb_points);
  }

  // loads grid file
  if (not(file_exists(filename_grid)) ||
      not(file_exists(filename_evaluated_grid))) {
    std::cout << "Grid file does not exist, has to be computed." << std::endl;
    intp = std::make_shared<IntpScalarUniformCubicSpline>(func, x1, x2,
                                                          error_bound);
    Vector_t grid{intp->get_grid_ref()};
    Vector_t evaluated_grid{intp->get_evaluated_grid_ref()};
    write_binary(filename_grid, grid);
    write_binary(filename_evaluated_grid, evaluated_grid);
  } else {
    std::cout << "Grid file exists and is read." << std::endl;
    Vector_t grid;
    Vector_t evaluated_grid;
    read_binary(filename_grid, grid);
    read_binary(filename_evaluated_grid, evaluated_grid);
    intp = std::make_shared<IntpScalarUniformCubicSpline>(grid, evaluated_grid);
  }
  std::cout
      << "interpolation of scalar radial contribution: interpolator grid size "
      << intp->get_grid_size() << std::endl;

  std::chrono::duration<double> elapsed{};
  auto start = std::chrono::high_resolution_clock::now();
  for (int j{0}; j < N_REPETITIONS; j++) {
    for (size_t i{0}; i < nb_iterations; i++) {
      points_tmp(i % nb_points) = intp->interpolate(points(i % nb_points));
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Interpolation of " << nb_iterations << " points"
            << " elapsed: " << elapsed.count() / N_REPETITIONS << " seconds"
            << std::endl;

  return 0;
}
