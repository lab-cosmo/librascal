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
      int max_radial{3};
      int max_angular{max_radial-1};
      json fc_hypers{
           {"type", "Constant"},
           {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
          };
      json hypers{{"gaussian_density", fc_hypers},
                {"max_radial", max_radial},
                {"max_angular", max_angular},
                {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
        };
      auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
      Vector_t tmp = Vector_t::Zero(fix.ref_points.size());
      for (auto _ : state) {
        for (size_t i{0}; i<fix.nb_iterations;i++) {
          benchmark::DoNotOptimize ( radial_contr.compute_contribution<AtomicSmearingType::Constant>(fix.ref_points(i % fix.nb_ref_points), 0.5) ) ;
        }
      }

      state.SetComplexityN(fix.nb_iterations);
      state.counters.insert({
          {"x1",fix.x1},
          {"x2",fix.x2},
          {"nb_iterations",fix.nb_iterations}
        });
    }

  template <class Fix>
  void BM_IntpRadCon(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    //std::cout <<  fix.log_error_bound << std::endl;
    int max_radial{3};
    int max_angular{max_radial-1};
    json fc_hypers{
         {"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
        };
    json hypers{{"gaussian_density", fc_hypers},
              {"max_radial", max_radial},
              {"max_angular", max_angular},
              {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
    };
    auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
    auto intp = InterpolatorVectorized<
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >();
    std::function<Matrix_t(double)> func = [&radial_contr](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x, 0.5);};
    intp.initialize(func, fix.x1, fix.x2, fix.error_bound); 
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_iterations;i++) {
        benchmark::DoNotOptimize( intp.interpolate(fix.ref_points(i % fix.ref_points.size())) );
      }
    }
    state.SetComplexityN(fix.nb_iterations);
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(error_bound)", fix.log_error_bound},
        {"log(max_error)", std::log10(intp.max_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",intp.grid.size()},
        {"random",fix.random}
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
        {"log(mean_error)",std::log10(fix.intp.mean_error)},
        {"log(max_error)",std::log10(fix.intp.max_error)},
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
        {"log(mean_error)",std::log10(fix.intp.mean_error)},
        {"log(max_error)",std::log10(fix.intp.max_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",fix.intp.grid.size()}
      });
  }

  // TODO(alex) I want to SetUp the function only once, but 
  // static object

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
        {"log(mean_error)",std::log10(fix.intp.mean_error)},
        {"log(max_error)",std::log10(fix.intp.max_error)},
        {"nb_iterations",fix.nb_iterations},
        {"grid_size",fix.intp.grid.size()}
      });
  }  
  
  //TwoDatasets one which you use for complexity measurements, thus all data is already initialized
  

  //auto intpfix_i_s{IntpFix<I_S>()};
  //BENCHMARK_CAPTURE(BM_IntpSimpFun, test_name, intpfix_i_s)->Apply(AllCombinationsArguments<I_S>)->Complexity();
  //BENCHMARK_TEMPLATE(BM_IntpSimpFun, IntpFix<I_B>)->Apply(AllCombinationsArguments<I_B>);
  
  // TODO(alex) since they are anyway global we can use BENCHMARK
  auto intp_fix{InterpolatorFixture<CFI_B>()};
  //BENCHMARK_CAPTURE(BM_RadCon, t, intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();
  //BENCHMARK_CAPTURE(BM_IntpRadCon, , intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();

  BENCHMARK_CAPTURE(BM_Hyp1f1, , intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();
  BENCHMARK_CAPTURE(BM_IntpHyp1f1, , intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();
  BENCHMARK_CAPTURE(BM_Hyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();
  BENCHMARK_CAPTURE(BM_IntpHyp1f1_DoNotOptimize, , intp_fix)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();



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
