/**
 * file  profile_vector_cubic_spline.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief used to profile the vector cubic spline method by storing the
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

#include <chrono>

#include "json.hpp"
#include "json_io.hh"
#include "profile_utility.hh"
#include "math/interpolator.hh"
#include "representations/representation_manager_spherical_expansion.hh"

using namespace rascal::math;  // NOLINT
using namespace rascal::internal;  // NOLINT

static constexpr int N_REPETITIONS = 200;
static constexpr int SEED = 1597463007;

/* Please execute this file one time before using the profiler to create the
 * grid.
 */
int main() {
  // radial contribution parameters
  int max_radial{3};
  int max_angular{max_radial - 1};
  json fc_hypers{{"type", "Constant"},
                 {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}};
  json hypers{
      {"gaussian_density", fc_hypers},
      {"max_radial", max_radial},
      {"max_angular", max_angular},
      {"compute_gradients", true},
      {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};
  auto radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
  std::function<Matrix_t(double)> func = [&radial_contr](double x) {
    return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x,
                                                                           0.5);
  };

  // interpolator parameters
  auto intp{InterpolatorVectorized<
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>>()};
  double x1{0};
  double x2{8};
  double mean_error_bound{1e-10};
  size_t nb_points = 1e6;
  size_t nb_iterations = 100000;
  const char * filename{"profile_vector_cubic_spline_grid.dat"};

  // loads grid file
  if (not(file_exists(filename))) {
    std::cout << "Grid file does not exist, has to be computed." << std::endl;
    intp.initialize(func, x1, x2, mean_error_bound);
    write_binary(filename, intp.grid);
  } else {
    std::cout << "Grid file exists and is read." << std::endl;
    Vector_t grid;
    read_binary(filename, grid);
    intp.initialize(func, grid);
  }
  std::cout << "vector interpolation of radial contribution: interpolator grid size " << intp.grid.size() << std::endl;

  // shuffle points to test interpolation method for uncorrelated requests
  srand(SEED); 
  Vector_t points_tmp = Vector_t::LinSpaced(nb_points, x1, x2);
  Vector_t points = Vector_t::Zero(nb_points);
  for (size_t i{0}; i < nb_points; i++) {
    points(i) = points_tmp(rand() % nb_points);
  }

  Matrix_t mat_tmp = Matrix_t::Zero(max_radial, max_angular + 1);
  std::chrono::duration<double> elapsed{};
  auto start = std::chrono::high_resolution_clock::now();
  for (int j{0}; j < N_REPETITIONS; j++) {
    for (size_t i{0}; i < nb_iterations; i++) {
      mat_tmp = intp.interpolate(points(i % nb_points));
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Interpolation of " << nb_iterations << " points"
            << " elapsed: " << elapsed.count() / N_REPETITIONS << " seconds"
            << std::endl;

  return 0;
}
