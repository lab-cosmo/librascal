/**
 * file  chrono_benchmark_interpolator.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief Benchmarks the interpolator using chrono, main usage for comparisment
 * with google benchmarks
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

#include "json_io.hh"
#include "profile_utility.hh"
#include "math/interpolator.hh"
#include "representations/representation_manager_spherical_expansion.hh"

using namespace rascal;  // NOLINT
using namespace rascal::math;  // NOLINT
using namespace rascal::internal;  // NOLINT

/* The number of repetitions to repeat the interpolation method for a number of
 * iterations. For the initialization process we only use one repititon, because
 * it takes comparatively long.
 */
static constexpr int N_REPETITIONS = 100;
static constexpr int SEED = 1597463007;

int main() {

    // Benchmark parameters //
    std::chrono::duration<double> elapsed{};
    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();


    // number of points of the grid used for profiling
    size_t nb_points = 1e6;
    // different numbers of iteration to estimate the scaling
    std::vector<size_t> nbs_iterations = {1000, 10000, 100000};

    // Interpolator parameters //
    auto intp{rascal::math::Interpolator<
        math::InterpolationMethod<math::InterpolationMethod_t::CubicSpline>,
        math::GridRational<math::GridType_t::Uniform, math::RefinementMethod_t::Exponential>,
        math::SearchMethod<math::SearchMethod_t::Uniform>>()};
    double x1{0};
    double x2{8};
    double error_bound{1e-5};
    std::function<double(double)> func;

    // shuffle points to test interpolation method for uncorrelated requests
    srand(SEED); 
    Vector_t points_tmp = Vector_t::LinSpaced(nb_points, x1, x2);
    Vector_t points = Vector_t::Zero(nb_points);
    for (size_t i{0}; i < nb_points; i++) {
      points(i) = points_tmp(rand() % nb_points);
    }

    ///* benchmark for the hyp1f1 function
    // *
    // * this benchmark compares the hyp1f1 function and its interpolation with
    // * the scalar interpolator. 
    // */

    //std::cout << "Start benchmark for hyp1f1" << std::endl;
    //// hyp1f1 function parameters //
    //double n = 5;
    //double l = 5;
    //double a = 0.5 * (n + l + 3);
    //double b = l + 1.5;
    //auto hyp1f1 = Hyp1f1(a, b, 200, 1e-15);
    //func = [&hyp1f1](double x) {
    //  return hyp1f1.calc(x);
    //};

    //// hyp1f1 //
    //for (size_t nb_iterations : nbs_iterations) {
    //  start = std::chrono::high_resolution_clock::now();
    //  for (int j{0}; j < N_REPETITIONS; j++) {
    //    for (size_t i{0}; i < nb_iterations; i++) {
    //      points_tmp(i % nb_points) = hyp1f1.calc(points(i % points.size()));
    //    }
    //  }
    //  finish = std::chrono::high_resolution_clock::now();

    //  elapsed = finish - start;
    //  std::cout << std::fixed;
    //  std::cout << "hyp1f1 elapsed: "
    //            << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                       .count() /
    //                   double(N_REPETITIONS)
    //            << " ns "
    //            << "hyp1f1 for " << nb_iterations << " points" << std::endl;
    //}
    //std::cout << std::endl;
    //
    //// interpolation of hyp1f1 //
    //
    //// measuring initialization process
    //start = std::chrono::high_resolution_clock::now();
    //intp.initialize(func, x1, x2, error_bound);
    //finish = std::chrono::high_resolution_clock::now();
    //elapsed = finish - start;
    //std::cout << std::fixed;
    //std::cout << "interpolation of hyp1f1: initialization elapsed time: "
    //          << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                     .count() 
    //          << " ns."
    //          << std::endl;
    //std::cout << "interpolation of hyp1f1: interpolator grid size " << intp.grid.size() << std::endl;
    //std::cout << std::endl;

    //// measuring interpolation method
    //for (size_t nb_iterations : nbs_iterations) {
    //  auto start = std::chrono::high_resolution_clock::now();
    //  for (int j{0}; j < N_REPETITIONS; j++) {
    //    for (size_t i{0}; i < nb_iterations; i++) {
    //      points_tmp(i % nb_points) = intp.interpolate(points(i % nb_points));
    //    }
    //  }
    //  auto finish = std::chrono::high_resolution_clock::now();

    //  elapsed = finish - start;
    //  std::cout << "interpolation of hyp1f1: interpolation elapsed time: "
    //            << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                       .count() /
    //                   double(N_REPETITIONS)
    //            << " ns "
    //            << "interpolation for " << nb_iterations << " points" << std::endl;
    //}

    //std::cout << std::endl << std::endl << std::endl;

    ///* benchmark for the scalar radial contribution
    // *
    // * this benchmark uses parameters for the radial contribution which
    // * make the radial contribution to a scalar function ℝ->ℝ by returning
    // * a 1x1 dimensional matrix. It is used for comparisment with the scalar
    // * interpolator.
    // */
    //
    //std::cout << "Start benchmark for scalar radial contribution" << std::endl;
    //// scalar radial contribution parameters
    int max_radial{1};
    int max_angular{0};
    json fc_hypers{{"type", "Constant"},
                   {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}};
    json hypers{
        {"gaussian_density", fc_hypers},
        {"max_radial", max_radial},
        {"max_angular", max_angular},
        {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};
    auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
    func = [&radial_contr](double x) {
      radial_contr.compute_neighbour_contribution(x, 0.5);
      return radial_contr.radial_integral_neighbour(0, 0);
    };

    //// scalar radial contribution //
    //for (size_t nb_iterations : nbs_iterations) {
    //  auto start = std::chrono::high_resolution_clock::now();
    //  for (int j{0}; j < N_REPETITIONS; j++) {
    //    for (size_t i{0}; i < nb_iterations; i++) {
    //      points_tmp(i % points.size()) =
    //          radial_contr.compute_contribution<AtomicSmearingType::Constant>(
    //              points(i % points.size()), 0.5)(0, 0);
    //    }
    //  }
    //  finish = std::chrono::high_resolution_clock::now();
    //  elapsed = finish - start;
    //  std::cout << std::fixed;
    //  std::cout << "scalar radial contribution elapsed: "
    //            << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                       .count() /
    //                   double(N_REPETITIONS)
    //            << " ns "
    //            << "for " << nb_iterations << " points" << std::endl;
    //}
    //std::cout << std::endl;


    //// interpolation of scalar radial contribution //
    //
    //// initialization
    //start = std::chrono::high_resolution_clock::now();
    //intp.initialize(func, x1, x2, error_bound);
    //finish = std::chrono::high_resolution_clock::now();
    //elapsed = finish - start;
    //std::cout << std::fixed;
    //std::cout << "interpolation of scalar radial contribution: initialization elapsed time: "
    //          << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                     .count()
    //          << " ns."
    //          << std::endl;
    //std::cout << "interpolation of scalar radial contribution: interpolator grid size " << intp.grid.size() << std::endl;
    //std::cout << std::endl;

    //// interpolation method
    //for (size_t nb_iterations : nbs_iterations) {
    //  start = std::chrono::high_resolution_clock::now();
    //  for (int j{0}; j < N_REPETITIONS; j++) {
    //    for (size_t i{0}; i < nb_iterations; i++) {
    //      points_tmp(i % nb_points) = intp.interpolate(points(i % nb_points));
    //    }
    //  }
    //  finish = std::chrono::high_resolution_clock::now();

    //  elapsed = finish - start;
    //  std::cout << "interpolation of scalar radial contribution: interpolation elapsed time: "
    //            << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                       .count() /
    //                   double(N_REPETITIONS)
    //            << " ns "
    //            << "for " << nb_iterations << " points" << std::endl;
    //}
    //std::cout << std::endl << std::endl << std::endl;

    ///* benchmark for the radial contribution
    // *
    // * This benchmark uses more realistic values for max_radial and max_angular
    // * resulting in a function which outputs a matrix and not a scalar version.
    // * It is used to for comparisment with the vector interpolator
    // *
    // */

    //std::cout << "Start benchmark for radial contribution" << std::endl;
    //// radial contribution parameters
    max_radial = 5;
    max_angular = max_radial;
    hypers = {{"gaussian_density", fc_hypers},
              {"max_radial", max_radial},
              {"max_angular", max_angular},
              {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};
    radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
    std::function<Matrix_t(double)> func_vec = [&radial_contr](double x) {
      return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x,
                                                                             0.5);
    };

    //// vector interpolator parameters 
    auto intp_vec{InterpolatorVectorized<
        InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
        GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
        SearchMethod<SearchMethod_t::Uniform>>()};

    //// radial contribution //
    //for (size_t nb_iterations : nbs_iterations) {
    //  start = std::chrono::high_resolution_clock::now();
    //  for (int j{0}; j < N_REPETITIONS; j++) {
    //    for (size_t i{0}; i < nb_iterations; i++) {
    //      Matrix_t points_vec_tmp =
    //          radial_contr.compute_contribution<AtomicSmearingType::Constant>(
    //              points(i % points.size()), 0.5);
    //    }
    //  }
    //  finish = std::chrono::high_resolution_clock::now();
    //  elapsed = finish - start;
    //  std::cout << std::fixed;
    //  std::cout << "radial contribution: elapsed time: "
    //            << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
    //                       .count() /
    //                   double(N_REPETITIONS)
    //            << " ns "
    //            << "for " << nb_iterations << " points" << std::endl;
    //}
    //std::cout << std::endl;

    // interpolation of radial contribution //

    // initialization
    start = std::chrono::high_resolution_clock::now();
    intp_vec.initialize(func, x1, x2, error_bound);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << std::fixed;
    std::cout << "vector interpolation of radial contribution: initialization elapsed time: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
                         .count()
              << " ns."
              << std::endl;
    std::cout << "vector interpolation of radial contribution: interpolator grid size " << intp_vec.grid.size() << std::endl;
    std::cout << std::endl;
    
    // TODO(alex) bug here do with debugger
    // iintpnterpolation method
    for (size_t nb_iterations : nbs_iterations) {
      start = std::chrono::high_resolution_clock::now();
      for (int j{0}; j < N_REPETITIONS; j++) {
        for (size_t i{0}; i < nb_iterations; i++) {
          Matrix_t points_vec_tmp =
              intp_vec.interpolate(points(i % nb_points));
        }
      }
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      std::cout << "vector interpolation of radial contribution: interpolation elapsed time: "
                << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
                           .count() /
                       double(N_REPETITIONS)
                << " ns "
                << "for " << nb_iterations << " points" << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;

    /* benchmark overhead
     *
     * measures the overhead not related to the interpolation method
     */
    std::cout << "Start benchmark for overhead" << std::endl;
    
    // measures the overhead of the measurement procedure //
    for (size_t nb_iterations : nbs_iterations) {
      start = std::chrono::high_resolution_clock::now();
      for (int j{0}; j < N_REPETITIONS; j++) {
        for (size_t i{0}; i < nb_iterations; i++) {
          // ~x65 faster than cubic spline
          Matrix_t points_vec_tmp = points(i % nb_points)
            * Matrix_t::Ones(max_radial, max_angular+1);
        }
      }
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      std::cout << "vector interpolation of radial contribution: \'measuring overhead\' elapsed time: "
                << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
                           .count() /
                       double(N_REPETITIONS)
                << " ns "
                << "" << nb_iterations << " iterations" << std::endl;
    }
    
    // measures overhead of the measurement procedure and the Interpolator //
    // class by using a method which does not only return an Identity*x //
    // matrix //
    for (size_t nb_iterations : nbs_iterations) {
      start = std::chrono::high_resolution_clock::now();
      for (int j{0}; j < N_REPETITIONS; j++) {
        for (size_t i{0}; i < nb_iterations; i++) {
           // ~x41 faster than cubic spline
           Matrix_t points_vec_tmp = intp_vec.interpolate_optimal(points(i %
             nb_points));
        }
      }
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      std::cout << "vector interpolation of radial contribution: \'measuring and interpolator class overhead\' elapsed time: "
                << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed)
                           .count() /
                       double(N_REPETITIONS)
                << " ns "
                << "" << nb_iterations << " iterations" << std::endl;
    }
    
    return 0;
}
