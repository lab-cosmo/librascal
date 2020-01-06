/**
 * @file  performance/profiles/profile_matrix_cubic_spline.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief used to profile the matrix cubic spline method by storing the
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

#include "rascal/math/interpolator.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/utils/json_io.hh"

#include <chrono>

static constexpr int N_REPETITIONS = 200;
static unsigned int SEED = 1597463007;

using namespace rascal;  // NOLINT

using math::InterpolatorMatrixUniformCubicSpline;
using math::Matrix_t;
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
  // radial contribution parameters
  int max_radial{3};
  int max_angular{max_radial};
  json fc_hypers{{"type", "Constant"},
                 {"gaussian_sigma", {{"value", 0.5}, {"unit", "AA"}}}};
  json hypers{
      {"gaussian_density", fc_hypers},
      {"max_radial", max_radial},
      {"max_angular", max_angular},
      {"compute_gradients", true},
      {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "AA"}}}}}};
  auto radial_contr =
      internal::RadialContribution<internal::RadialBasisType::GTO>(hypers);
  std::function<Matrix_t(double)> func = [&radial_contr](double x) {
    return radial_contr
        .compute_contribution<internal::AtomicSmearingType::Constant>(x, 0.5);
  };

  // interpolator parameters
  using IntpVectorUniformCubicSpline =
      InterpolatorMatrixUniformCubicSpline<RefinementMethod_t::Exponential>;
  std::shared_ptr<IntpVectorUniformCubicSpline> intp;
  double x1{0};
  double x2{8};
  double error_bound{1e-10};
  size_t nb_points = 1e6;
  size_t nb_iterations = 100000;
  const char * filename_grid{"profile_vector_cubic_spline_grid.grid"};
  const char * filename_evaluated_grid{
      "profile_vector_cubic_spline_grid.evaluated_grid"};

  // loads grid file
  if (not(rascal::file_exists(filename_grid)) ||
      not(rascal::file_exists(filename_evaluated_grid))) {
    std::cout << "Grid file does not exist, has to be computed." << std::endl;
    Matrix_t result = func(x1);
    int cols{static_cast<int>(result.cols())};
    int rows{static_cast<int>(result.rows())};
    intp = std::make_shared<IntpVectorUniformCubicSpline>(
        func, x1, x2, error_bound, cols, rows);
    Vector_t grid{intp->get_grid_ref()};
    Matrix_t evaluated_grid{intp->get_evaluated_grid_ref()};
    write_binary(filename_grid, grid);
    write_binary(filename_evaluated_grid, evaluated_grid);
  } else {
    std::cout << "Grid file exists and is read." << std::endl;
    Vector_t grid;
    Matrix_t evaluated_grid;
    read_binary(filename_grid, grid);
    read_binary(filename_evaluated_grid, evaluated_grid);
    std::cout << grid.size() << " " << evaluated_grid.cols() << " "
              << evaluated_grid.rows() << std::endl;
    intp = std::make_shared<IntpVectorUniformCubicSpline>(
        grid, evaluated_grid, max_radial, max_angular + 1);
  }
  std::cout
      << "vector interpolation of radial contribution: interpolator grid size "
      << intp->get_grid_size() << std::endl;

  // shuffle points to test interpolation method for uncorrelated requests
  srand(SEED);
  Vector_t points_tmp = Vector_t::LinSpaced(nb_points, x1, x2);
  Vector_t points = Vector_t::Zero(nb_points);
  for (size_t i{0}; i < nb_points; i++) {
    points(i) = points_tmp(rand_r(&SEED) % nb_points);
  }

  Matrix_t mat_tmp = Matrix_t::Zero(max_radial, max_angular + 1);
  std::chrono::duration<double> elapsed{};
  auto start = std::chrono::high_resolution_clock::now();
  for (int j{0}; j < N_REPETITIONS; j++) {
    for (size_t i{0}; i < nb_iterations; i++) {
      mat_tmp = intp->interpolate(points(i % nb_points));
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Interpolation of " << nb_iterations << " points"
            << " elapsed: " << elapsed.count() / N_REPETITIONS << " seconds"
            << std::endl;

  return 0;
}
