#include "benchmark_interpolator.hh"


using namespace rascal::math;

namespace rascal {
  namespace internal {
    ////////////
    // Issues //
    ////////////
    // Currently there is not way to print out the results in scientific notation which make it harder to interpret the results. There is a pull request takling this issue https://github.com/google/benchmark/pull/821 but it is not yet merged
    
    // TODO(alex) To have a credible time benchmarks for the initialization function.
    // distance function caluclates for (radial_max, angular_max+1) points, this must be incuded
    // in the intp, maybe make extra class MatrixIntp

    //test interpoltar matrix
    //1.tests Interpoltaro for gauss function compares only initialization speed


    template <class Fix>
    void BM_RadCon(benchmark::State& state, Fix & fix) {
      fix.SetUp(state);
      int max_radial{1};
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
      Vector_t tmp = Vector_t::Zero(fix.points.size());
      for (auto _ : state) {
        for (size_t i{0}; i<fix.nb_points;i++) {
          tmp(i % fix.points.size()) = radial_contr.compute_contribution<AtomicSmearingType::Constant>(fix.points(i % fix.points.size()), 0.5)(0,0);
        }
      }

      state.SetComplexityN(fix.nb_points);
      state.counters.insert({
          {"x1",fix.x1},
          {"x2",fix.x2},
          {"nb_points",fix.nb_points}
        });
    }

  template <class Fix>
  void BM_IntpRadCon(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    //std::cout <<  fix.log_mean_error_bound << std::endl;
    int max_radial{1};
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
    using AdaptiveInterpolator = Interpolator <
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;
    const int nb_intps{max_radial*(max_angular+1)};
    //std::vector<std::function<double(double)>> funcs(nb_intps); 
    std::vector<AdaptiveInterpolator> intps(nb_intps);
    for (int i{0}; i<nb_intps;i++) { 
      //std::cout<< "Iterating=" << i << std::endl;
      auto func = [&radial_contr, nb_intps, max_angular, max_radial,i](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x, 0.5)(i/(max_angular+1),i % max_radial);};
      intps.push_back(AdaptiveInterpolator());
      intps.at(i).initialize(func, fix.x1, fix.x2, fix.mean_error_bound); 
    }
    Vector_t tmp = Vector_t::Zero(fix.points.size());
    auto intp = intps.at(0);
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_points;i++) {
        tmp(i % fix.points.size()) = intp.interpolate(fix.points(i));
      }
    }
    //for (auto _ : state) {
    //  for (int j{0}; j<nb_intps; j++) {
    //    for (int i{0}; i<fix.points.size();i++) {
    //      intps.at(j).interpolate(fix.points(i));
    //    }
    //  }
    //}
    //fix.TearDown();
    state.SetComplexityN(fix.nb_points);
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"nb_points",fix.nb_points},
        {"grid_size",intps.at(0).grid.size()},
        {"random",fix.random}
      });
  }



  // Difference between writing data into an array and DoNotOptimize is 1% (Writing is 1% faster
  // benchmark::ClobberMemory(); does not change anything
  template <class Fix>
  void BM_Hyp1f1(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    double n = 10;
    double l = 10;
    double a = 0.5*(n+l+3);
    double b = l+1.5;
    auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
    Vector_t tmp = Vector_t::Zero(fix.points.size());
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_points;i++) {
        //benchmark::DoNotOptimize (
        //  hyp1f1.calc(fix.points(i % fix.points.size()))
        //);
        tmp(i % fix.points.size()) = hyp1f1.calc(fix.points(i % fix.points.size()));
        
      }
    }
    state.SetComplexityN(fix.nb_points);
    state.counters.insert({
        {"nb_points",fix.nb_points}
      });
  }

  template <class Fix>
  void BM_IntpHyp1f1(benchmark::State &state, Fix & fix) {
    fix.SetUp(state);
    auto intp{Interpolator <
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >()};
    double n = 10;
    double l = 10;
    double a = 0.5*(n+l+3);
    double b = l+1.5;
    auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
    auto func = [&hyp1f1](double x) {return hyp1f1.calc(x);};
    intp.initialize(func, fix.x1, fix.x2, fix.mean_error_bound); 
    Vector_t tmp = Vector_t::Zero(fix.points.size());
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_points;i++) {
        //benchmark::DoNotOptimize (
        //  intp.interpolate(fix.points(i % fix.points.size()))
        //);
        tmp(i % fix.points.size()) = intp.interpolate(fix.points(i % fix.points.size()));
      }
    }
    state.SetComplexityN(fix.nb_points);
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"log(mean_error)",std::log10(intp.mean_error)},
        //{"log(max_error)",std::log10(intp.max_error)},
        {"nb_points",fix.nb_points},
        {"grid_size",intp.grid.size()}
      });
  }

  // TODO(alex) I want to SetUp the function only once, but 
  // static object

  template <class IntpFix>
  void BM_IntpSimpFun(benchmark::State &state, IntpFix & fix) {
    fix.SetUp(state);
    //std::cout << "Start iteration " << fix.nb_points << std::endl;
    for (auto _ : state) {
      for (size_t i{0}; i<fix.nb_points;i++) {
        benchmark::DoNotOptimize (
          fix.intp.interpolate(fix.points(i % fix.points.size()))
        );
      }
    }
    //std::cout << "End iteration" << std::endl;
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"log(mean_error)",std::log10(fix.intp.mean_error)},
        {"log(max_error)",std::log10(fix.intp.max_error)},
        {"nb_points",fix.nb_points},
        {"grid_size",fix.intp.grid.size()}
      });
  }  
  
  //TwoDatasets one which you use for complexity measurements, thus all data is already initialized
  

  //auto intpfix_i_s{IntpFix<I_S>()};
  //BENCHMARK_CAPTURE(BM_IntpSimpFun, test_name, intpfix_i_s)->Apply(AllCombinationsArguments<I_S>)->Complexity();
  //BENCHMARK_TEMPLATE(BM_IntpSimpFun, IntpFix<I_B>)->Apply(AllCombinationsArguments<I_B>);
  
  // TODO(alex) since they are anyway global we can use BENCHMARK
  auto intpfix_cfi_b_r1{IntpFix<CFI_B>()};
  BENCHMARK_CAPTURE(BM_RadCon, t, intpfix_cfi_b_r1)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();
  auto intpfix_cfi_b_r2{IntpFix<CFI_B>()};
  BENCHMARK_CAPTURE(BM_IntpRadCon, t, intpfix_cfi_b_r2)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();

  auto intpfix_cfi_b_h2{IntpFix<CFI_B>()};
  BENCHMARK_CAPTURE(BM_Hyp1f1, t2, intpfix_cfi_b_h2)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();

  auto intpfix_cfi_b_h1{IntpFix<CFI_B>()};
  BENCHMARK_CAPTURE(BM_IntpHyp1f1, t, intpfix_cfi_b_h1)->Apply(AllCombinationsArguments<CFI_B>)->Complexity();



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
