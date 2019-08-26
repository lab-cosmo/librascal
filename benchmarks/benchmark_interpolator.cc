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
// TODO(alex) still needed?
#include "representations/feature_manager_block_sparse.hh"


using namespace rascal::math;

namespace rascal {
  namespace internal {
    template <class Fix>
    void BM_SphExp(benchmark::State& state, Fix & fix) {
      fix.SetUp(state);
      for (auto _ : state) {
        fix.representation_ptr->compute();
      }
      // TODO(alex) I would like to print interpolator information, but it is covered under a lot of layers
      state.counters.insert({
          {"max_radial", fix.max_radial},
          {"max_angular", fix.max_radial},
          {"cutoff", fix.cutoff},
          {"nb_neighbours", fix.nb_neighbours}
        });
    }

    template <class Fix>
    void BM_RadCon(benchmark::State& state, Fix & fix) {
      fix.SetUp(state);
      for (auto _ : state) {
        for (size_t i{0}; i<fix.nb_iterations;i++) {
          benchmark::DoNotOptimize( fix.func(fix.ref_points(i % fix.nb_ref_points)) );
        }
      }
      state.SetComplexityN(fix.nb_iterations);
      state.counters.insert({
          {"x1",fix.x1},
          {"x2",fix.x2},
          {"nb_iterations",fix.nb_iterations},
          {"max_radial", fix.max_radial}
        });
    }

    template <class Fix>
    void BM_IntpRadCon(benchmark::State &state, Fix & fix) {
      fix.SetUp(state);
      for (auto _ : state) {
        for (size_t i{0}; i<fix.nb_iterations;i++) {
          benchmark::DoNotOptimize( fix.intp.interpolate(fix.ref_points(i % fix.ref_points.size())) );
        }
      }
      state.SetComplexityN(fix.nb_iterations);
      state.counters.insert({
          {"x1",fix.x1},
          {"x2",fix.x2},
          {"nb_iterations",fix.nb_iterations},
          {"max_radial", fix.max_radial},
          {"log(error_bound)", fix.log_error_bound},
          {"log(max_grid_error)", std::log10(fix.intp.max_grid_error)},
          {"grid_size", fix.intp.grid.size()},
        });
    }



  // Difference between writing data into an array and DoNotOptimize is 1% (Writing is 1% faster
  template <class Fix>
  void BM_Hyp1f1(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        tmp(i % fix.ref_points.size()) = fix.func(fix.ref_points(i % fix.ref_points.size()));
        
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({
        {"nb_iterations",fix.nb_iterations}
      });
  }
  template <class Fix>
  void BM_Hyp1f1_DoNotOptimize(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        benchmark::DoNotOptimize (
          fix.func(fix.ref_points(i % fix.ref_points.size()))
        );
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({
        {"nb_iterations",fix.nb_iterations}
      });
  }

  template <class Fix>
  void BM_IntpHyp1f1(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    // to prevent optimization
    Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        tmp(i % fix.ref_points.size()) = fix.intp.interpolate(fix.ref_points(i % fix.ref_points.size()));
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(error_bound)", fix.log_error_bound},
        {"log(mean_grid_error)",std::log10(fix.intp.mean_grid_error)},
        {"log(max_grid_error)",std::log10(fix.intp.max_grid_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",fix.intp.grid.size()}
      });
  }
  template <class Fix>
  void BM_IntpHyp1f1_DoNotOptimize(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    // to prevent optimization
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        benchmark::DoNotOptimize (
          fix.intp.interpolate(fix.ref_points(i % fix.ref_points.size()))
        );
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(error_bound)", fix.log_error_bound},
        {"log(mean_grid_error)",std::log10(fix.intp.mean_grid_error)},
        {"log(max_grid_error)",std::log10(fix.intp.max_grid_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",fix.intp.grid.size()}
      });
  }


  template <class IntpFix>
  void BM_IntpSimpFun(benchmark::State &state, IntpFix & fix) {
    fix.SetUp(state);
    //std::cout << "Start iteration " << fix.nb_iterations << std::endl;
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        benchmark::DoNotOptimize (
          fix.intp.interpolate(fix.ref_points(i % fix.ref_points.size()))
        );
      }
    }
    //std::cout << "End iteration" << std::endl;
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(error_bound)", fix.log_error_bound},
        {"log(mean_grid_error)",std::log10(fix.intp.mean_grid_error)},
        {"log(max_grid_error)",std::log10(fix.intp.max_grid_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",fix.intp.grid.size()}
      });
  }  
  


  //auto sph_exp_intp_fix = SphericalExpansionBFixture<SphericalExpansionDataset>(true, false);
  //BENCHMARK_CAPTURE(BM_SphExp, use_intp_no_gradient, sph_exp_intp_fix)->Apply(AllCombinationsArguments<SphericalExpansionDataset>);
  //auto sph_exp_fix = SphericalExpansionBFixture<SphericalExpansionDataset>(false, false);
  //BENCHMARK_CAPTURE(BM_SphExp, no_intp_no_gradient, sph_exp_fix)->Apply(AllCombinationsArguments<SphericalExpansionDataset>);

  //auto sph_exp_intp_gradient_fix = SphericalExpansionBFixture<SphericalExpansionDataset>(true, true);
  //BENCHMARK_CAPTURE(BM_SphExp, use_intp_comp_gradient , sph_exp_intp_gradient_fix)->Apply(AllCombinationsArguments<SphericalExpansionDataset>);
  auto sph_exp_gradient_fix = SphericalExpansionBFixture<SphericalExpansionDataset>(false, true);
  BENCHMARK_CAPTURE(BM_SphExp, no_intp_comp_gradient, sph_exp_gradient_fix)->Apply(AllCombinationsArguments<SphericalExpansionDataset>);

  //auto intp_vec_fix{InterpolatorVectorizedFixture<RadConDataset>()};
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  
  //auto intp_vec_fix1{InterpolatorVectorizedFixture<RadConDataset1>()};
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix1)->Apply(AllCombinationsArguments<RadConDataset1>)->Complexity();
  //BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix1)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  //auto intp_vec_fix2{InterpolatorVectorizedFixture<RadConDataset2>()};
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix2)->Apply(AllCombinationsArguments<RadConDataset2>)->Complexity();
  //BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix2)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  //auto intp_vec_fix3{InterpolatorVectorizedFixture<RadConDataset3>()};
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix3)->Apply(AllCombinationsArguments<RadConDataset3>)->Complexity();
  //BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix3)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();

  //auto intp_fix{InterpolatorScalarFixture<Hyp1f1Dataset>()};
  //BENCHMARK_CAPTURE(BM_Hyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<Hyp1f1Dataset>)->Complexity();
  //BENCHMARK_CAPTURE(BM_IntpHyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<Hyp1f1Dataset>)->Complexity();
 
  } // namespace internal
} // namespace rascal

BENCHMARK_MAIN();
