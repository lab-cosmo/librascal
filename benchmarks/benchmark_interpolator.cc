/**
 * file  benchmark_interpolator.cc
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

  // Hyp1f1 function
  template <class BFixture>
  void BM_Hyp1f1(benchmark::State & state, BFixture & fix) {
    fix.SetUp(state);
    // to prevent optimization we copy the results
    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp(i % fix.ref_points.size()) =
            fix.func(fix.ref_points(i % fix.ref_points.size()));
      }
      /**
       * Alternatively the DoNotOptimize function can be used, however I am
       * not sure how much is not optimized.
       *
       * for (size_t i{0}; i < fix.nb_iterations; i++) {
       *  benchmark::DoNotOptimize(
       *      fix.func(fix.ref_points(i % fix.ref_points.size())));
       *}
       */
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({{"nb_iterations", fix.nb_iterations}});
  }

  // Hyp1f1 using the interpolator
  template <class BFixture>
  void BM_Hyp1f1Intp(benchmark::State & state, BFixture & fix) {
    fix.SetUp(state);
    // to prevent optimization
    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp(i % fix.ref_points.size()) =
            fix.intp->interpolate(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    json info = fix.intp->compute_interpolator_information();

    state.counters.insert({{"x1", fix.x1},
                           {"x2", fix.x2},
                           {"log(error_bound)", fix.log_error_bound},
                           {"log(mean_grid_error)",
                            std::log10(info["mean_grid_error"].get<double>())},
                           {"log(max_grid_error)",
                            std::log10(info["max_grid_error"].get<double>())},
                           {"nb_iterations", fix.nb_iterations},
                           {"grid_size", info["grid_size"]}});
  }

  template <class Fix>
  void BM_RadCon(benchmark::State & state, Fix & fix) {
    fix.SetUp(state);
    Matrix_t tmp = Matrix_t::Zero(fix.max_radial, fix.max_angular + 1);
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp = fix.func(fix.ref_points(i % fix.nb_ref_points));
      }
    }
  }

  // Hyp1f1 using the interpolator
  template <class BFixture>
  void BM_RadConIntp(benchmark::State & state, BFixture & fix) {
    fix.SetUp(state);
    // to prevent optimization
    Matrix_t tmp = Matrix_t::Zero(fix.max_radial, fix.max_angular + 1);
    for (auto _ : state) {
      for (size_t i{0}; i < fix.nb_iterations; i++) {
        tmp = fix.intp->interpolate(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    json info = fix.intp->compute_interpolator_information();

    state.counters.insert({{"nb_iterations", fix.nb_iterations},
                           {"max_radial", fix.max_radial},
                           {"max_angular", fix.max_angular},
                           {"x1", fix.x1},
                           {"x2", fix.x2},
                           {"log(error_bound)", fix.log_error_bound},
                           {"log(mean_grid_error)",
                            std::log10(info["mean_grid_error"].get<double>())},
                           {"log(max_grid_error)",
                            std::log10(info["max_grid_error"].get<double>())},
                           {"grid_size", info["grid_size"]}});
  }

  template <class BFixture>
  void BM_SphExp(benchmark::State & state, BFixture & fix) {
    fix.SetUp(state);
    for (auto _ : state) {
      fix.representation_ptr->compute(fix.manager);
    }
    // TODO(alex) I would like to print interpolator information, but it is
    // covered under a lot of layers, not sure if
    state.counters.insert({{"max_radial", fix.max_radial},
                           {"max_angular", fix.max_angular},
                           {"cutoff", fix.cutoff},
                           {"nb_neighbours", fix.nb_neighbours}});
  }

  /**
   * Hyp1f1 for the scalar interpolator
   */
  auto intp_fix{InterpolatorScalarBFixture<Hyp1f1Dataset>()};
  BENCHMARK_CAPTURE(BM_Hyp1f1, , intp_fix)
      ->Apply(AllCombinationsArguments<Hyp1f1Dataset>)
      ->Complexity();
  BENCHMARK_CAPTURE(BM_Hyp1f1Intp, , intp_fix)
      ->Apply(AllCombinationsArguments<Hyp1f1Dataset>)
      ->Complexity();

  /**
   * RadialContribution for the matrix interpolator
   */
  auto intp_mat_fix{InterpolatorMatrixBFixture<RadialContributionDataset>()};
  BENCHMARK_CAPTURE(BM_RadCon, , intp_mat_fix)
      ->Apply(AllCombinationsArguments<RadialContributionDataset>)
      ->Complexity();
  BENCHMARK_CAPTURE(BM_RadConIntp, , intp_mat_fix)
      ->Apply(AllCombinationsArguments<RadialContributionDataset>)
      ->Complexity();

  /**
   * Spherical Expansion without gradient benchmarks
   */
  auto sph_exp_fix =
      SphericalExpansionBFixture<SphericalExpansionDataset>(false, false);
  BENCHMARK_CAPTURE(BM_SphExp, no_intp_no_gradient, sph_exp_fix)
      ->Apply(AllCombinationsArguments<SphericalExpansionDataset>)
      ->Complexity();
  auto sph_exp_intp_fix =
      SphericalExpansionBFixture<SphericalExpansionDataset>(true, false);
  BENCHMARK_CAPTURE(BM_SphExp, use_intp_no_gradient, sph_exp_intp_fix)
      ->Apply(AllCombinationsArguments<SphericalExpansionDataset>)
      ->Complexity();

  /**
   * Spherical Expansion with gradient benchmarks
   */
  auto sph_exp_gradient_fix =
      SphericalExpansionBFixture<SphericalExpansionDataset>(false, true);
  BENCHMARK_CAPTURE(BM_SphExp, no_intp_comp_gradient, sph_exp_gradient_fix)
      ->Apply(AllCombinationsArguments<SphericalExpansionDataset>)
      ->Complexity();
  auto sph_exp_intp_gradient_fix =
      SphericalExpansionBFixture<SphericalExpansionDataset>(true, true);
  BENCHMARK_CAPTURE(BM_SphExp, use_intp_comp_gradient,
                    sph_exp_intp_gradient_fix)
      ->Apply(AllCombinationsArguments<SphericalExpansionDataset>)
      ->Complexity();

}  // namespace rascal

BENCHMARK_MAIN();
