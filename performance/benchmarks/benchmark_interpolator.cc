/**
 * @file  benchmark_interpolator.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief benchmarks for the interpolator
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

#include "benchmark_interpolator.hh"

namespace rascal {
  /**
   * To compare the benchmark results with the chrono library, these benchmarks
   * have been also created with chrono and can be found in
   * https://gist.github.com/agoscinski/7d0e812041bea6b1cfcd5b44061d9a9f
   */

  // Benchmark for Hyp1f1 without interpolator
  template <class BFixture>
  void bm_hyp1f1(benchmark::State & state, BFixture & fix) {
    fix.setup(state);

    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      /**
       * To prevent that this step is optimized away from the compiler we copy
       * the results to some temporary, the DoNotOptimize function can be used
       * alternatively, however I am not sure how much is not optimized. One
       * gets results deviating by around 10%.
       *
       * for (size_t i{0}; i < fix.nb_iterations; i++) {
       *  benchmark::DoNotOptimize(
       *      fix.func(fix.ref_points(i % fix.ref_points.size())));
       * }
       */
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp(i % fix.ref_points.size()) =
            fix.func(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({{"nb_iterations", fix.nb_iterations}});
  }

  // Benchmark for Hyp1f1 with interpolator
  template <class BFixture>
  void bm_hyp1f1_intp(benchmark::State & state, BFixture & fix) {
    fix.setup(state);
    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        // To prevent that this step is optimized away from the compiler we copy
        // the results to some temporary.
        tmp(i % fix.ref_points.size()) =
            fix.intp->interpolate(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    Vector_t pointwise_grid_error =
        compute_pointwise_absolute_grid_error(fix.intp);
    double mean_grid_error{pointwise_grid_error.mean()};
    double max_grid_error{pointwise_grid_error.maxCoeff()};
    int grid_size{fix.intp->get_grid_size()};

    state.counters.insert(
        {{"x1", fix.x1},
         {"x2", fix.x2},
         {"log(error_bound)", fix.log_error_bound},
         {"log(mean_grid_error)", std::log10(mean_grid_error)},
         {"log(max_grid_error)", std::log10(max_grid_error)},
         {"nb_iterations", fix.nb_iterations},
         {"grid_size", grid_size}});
  }

  // Benchmark for RadialContribution without interpolator
  template <class Fix>
  void bm_radial_contr(benchmark::State & state, Fix & fix) {
    fix.setup(state);
    Matrix_t tmp = Matrix_t::Zero(fix.max_radial, fix.max_angular + 1);
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp = fix.func(fix.ref_points(i % fix.nb_ref_points));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
  }

  // Benchmark for RadialContribution with interpolator
  template <class BFixture>
  void bm_radial_contr_intp(benchmark::State & state, BFixture & fix) {
    fix.setup(state);
    // to prevent optimization
    Matrix_t tmp = Matrix_t::Zero(fix.max_radial, fix.max_angular + 1);
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp = fix.intp->interpolate(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    Matrix_t pointwise_grid_error =
        compute_pointwise_absolute_grid_error(fix.intp);
    double mean_grid_error{pointwise_grid_error.mean()};
    double max_grid_error{pointwise_grid_error.maxCoeff()};
    int grid_size{fix.intp->get_grid_size()};

    state.counters.insert(
        {{"nb_iterations", fix.nb_iterations},
         {"max_radial", fix.max_radial},
         {"max_angular", fix.max_angular},
         {"x1", fix.x1},
         {"x2", fix.x2},
         {"log(error_bound)", fix.log_error_bound},
         {"log(mean_grid_error)", std::log10(mean_grid_error)},
         {"log(max_grid_error)", std::log10(max_grid_error)},
         {"grid_size", grid_size}});
  }

  // Benchmark for SphericalExpansion with or without the interpolator
  template <class BFixture>
  void bm_spherical(benchmark::State & state, BFixture & fix) {
    fix.setup(state);
    for (auto _ : state) {
      fix.representation_ptr->compute(fix.manager);
    }
    // TODO(alex) I would like to print interpolator information, but it is
    // covered under a lot of layers, not sure if
    state.SetComplexityN(fix.max_radial * (fix.max_angular + 1));
    state.counters.insert({{"max_radial", fix.max_radial},
                           {"max_angular", fix.max_angular},
                           {"cutoff", fix.cutoff},
                           {"accuracy", std::log10(fix.error_bound)},
                           {"nb_neighbours", fix.nb_neighbours}});
  }

  /**
   * Hyp1f1 for the scalar interpolator
   *
   * auto intp_fix{InterpolatorScalarBFixture<Hyp1f1Dataset>()};
   * BENCHMARK_CAPTURE(bm_hyp1f1, , intp_fix)
   *     ->Apply(all_combinations_of_arguments<Hyp1f1Dataset>)
   *     ->Complexity();
   * BENCHMARK_CAPTURE(bm_hyp1f1_intp, , intp_fix)
   *     ->Apply(all_combinations_of_arguments<Hyp1f1Dataset>)
   *     ->Complexity();
   */

  /**
   * RadialContribution for the matrix interpolator
   */
  auto intp_mat_fix{InterpolatorMatrixBFixture<RadialContributionDataset>()};
  BENCHMARK_CAPTURE(bm_radial_contr, , intp_mat_fix)
      ->Apply(all_combinations_of_arguments<RadialContributionDataset>)
      ->Complexity();
  BENCHMARK_CAPTURE(bm_radial_contr_intp, , intp_mat_fix)
      ->Apply(all_combinations_of_arguments<RadialContributionDataset>)
      ->Complexity();

  /**
   * Spherical Expansion without gradient benchmarks
   */
  auto sph_expansion_fix =
      SphericalExpansionBFixture<SphericalDataset>(false, false);
  BENCHMARK_CAPTURE(bm_spherical, expansion_no_intp_no_gradient,
                    sph_expansion_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();
  auto sph_expansion_intp_fix =
      SphericalExpansionBFixture<SphericalDataset>(true, false);
  BENCHMARK_CAPTURE(bm_spherical, expansion_use_intp_no_gradient,
                    sph_expansion_intp_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();

  /**
   * Spherical Expansion with gradient benchmarks
   */
  auto sph_expansion_gradient_fix =
      SphericalExpansionBFixture<SphericalDataset>(false, true);
  BENCHMARK_CAPTURE(bm_spherical, expansion_no_intp_comp_gradient,
                    sph_expansion_gradient_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();
  auto sph_expansion_intp_gradient_fix =
      SphericalExpansionBFixture<SphericalDataset>(true, true);
  BENCHMARK_CAPTURE(bm_spherical, expansion_use_intp_comp_gradient,
                    sph_expansion_intp_gradient_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();

  /**
   * Spherical Invariants without gradient benchmarks
   */
  auto sph_inv_fix =
      SphericalInvariantsBFixture<SphericalDataset>(false, false);
  BENCHMARK_CAPTURE(bm_spherical, invariant_no_intp_no_gradient, sph_inv_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();
  auto sph_inv_intp_fix =
      SphericalInvariantsBFixture<SphericalDataset>(true, false);
  BENCHMARK_CAPTURE(bm_spherical, invariant_use_intp_no_gradient,
                    sph_inv_intp_fix)
      ->Apply(all_combinations_of_arguments<SphericalDataset>)
      ->Complexity();

}  // namespace rascal

BENCHMARK_MAIN();
