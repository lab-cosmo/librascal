#include "benchmark_interpolator.hh"


using namespace rascal::math;

namespace rascal {
  namespace internal {
    // 


    ////////////
    // Issues //
    ////////////
    // TODO(alex) To have a credible time benchmarks for the initialization function.
    // distance function caluclates for (radial_max, angular_max+1) points, this must be incuded
    // in the intp, maybe make extra class MatrixIntp

    //test interpoltar matrix
    //1.tests Interpoltaro for gauss function compares only initialization speed


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
  
  
  //auto intp_vec_fix{InterpolatorVectorizedFixture<RadConDataset>()};
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  auto intp_vec_fix1{InterpolatorVectorizedFixture<RadConDataset1>()};
  BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix1)->Apply(AllCombinationsArguments<RadConDataset1>)->Complexity();
  BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix1)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  auto intp_vec_fix2{InterpolatorVectorizedFixture<RadConDataset2>()};
  BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix2)->Apply(AllCombinationsArguments<RadConDataset2>)->Complexity();
  BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix2)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();
  auto intp_vec_fix3{InterpolatorVectorizedFixture<RadConDataset3>()};
  BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_vec_fix3)->Apply(AllCombinationsArguments<RadConDataset3>)->Complexity();
  BENCHMARK_CAPTURE(BM_RadCon, , intp_vec_fix3)->Apply(AllCombinationsArguments<RadConDataset>)->Complexity();

  auto intp_fix{InterpolatorFixture<Hyp1f1Dataset>()};
  BENCHMARK_CAPTURE(BM_Hyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<Hyp1f1Dataset>)->Complexity();
  BENCHMARK_CAPTURE(BM_IntpHyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<Hyp1f1Dataset>)->Complexity();



  ///////////////////////////////////////////////
  // PSEUDO CODE HOW TO USE TESTS FOR MANAGERS //
  ///////////////////////////////////////////////

  //class ManagerDataFixture {
  // Manager_t = Strict;
  // public:
  //  
  //  static const json data() {
  //    return {
  //      {"filenames", {"crystal.json", "structure.json"}},
  //      {"sigmas", {0.1,0.3,0.5}}
  //      };
  //  }
  //}

  //template <class ManagerDataFixture>
  //class ManagerFixture : public ::benchmark::Fixture {
  // public:
  //  using Manager_t = DataFixture::Manager_t;
  //  // This is executed one time when this benchmark is used
  //  void SetUp(const ::benchmark::State& state) {
  //    json data = ManagerDataFixture::data();
  //    for (json::iterator it = data.begin(); it != data.end(); ++it) {
  //      std::string filename = data["filenames"].at(i).get<string>();
  //      auto manager = create_manager<Manager_t>(filename);
  //      managers.push_back();        
  //    }
  //  }
  //  static const json data() {
  //    return {
  //      {"filenames", {"crystal.json", "structure.json"}},
  //      {"max_angular", {5,10,15}}
  //      };
  //  }
  //  Manager_t managers;
  //  std::vector<size_t> max_angulars;
  //}

  //BENCHMARK_DEFINE_F(ManagerFixture<ManagerDataFixture>, ManagerBench)(benchmark::State &state) {
  //  // This is executed #repetitions
  //  auto manager = managers.at(state.range(0)); 
  //  size_t max_angular = max_angulars.at(state.range(1)); 
  //  for (auto _ : state) {
  //    // This is executed #iterations
  //  }
  //}
  //BENCHMARK_DEFINE_F(ManagerFixture<ManagerDataFixture>, ManagerBench)->Apply(AllCombinationsArguments<ManagerDataFixture>)


 //https://github.com/google/benchmark/issues/387
 
  } // namespace internal
} // namespace rascal

BENCHMARK_MAIN();
