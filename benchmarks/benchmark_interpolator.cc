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
    void BM_RadialContr(benchmark::State& state) {
      auto fix = Fix(state);
      int max_radial{2};
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
      for (auto _ : state) {
        for (int i{0}; i<fix.points.size();i++) {
          radial_contr.compute_contribution<AtomicSmearingType::Constant>(fix.points(i), 0.5);
        }
      }
      state.counters.insert({
          {"x1",fix.x1},
          {"x2",fix.x2},
          {"nb_points",fix.nb_points}
        });
    }

  template <class Fix>
  void BM_InterpolatorRadialContr(benchmark::State &state) {
    auto fix = Fix(state);
    int max_radial{2};
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
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
      SearchMethod<SearchMethod_t::Hunt>
        >;
    const int nb_intps{max_radial*(max_angular+1)};
    //std::vector<std::function<double(double)>> funcs(nb_intps); 
    std::vector<AdaptiveInterpolator> intps(nb_intps);
    for (int i{0}; i<nb_intps;i++) { 
      //std::cout<< "Iterating=" << i << std::endl;
      auto func = [&radial_contr, nb_intps, max_angular, max_radial,i](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x, 0.5)(i/(max_angular+1),i % max_radial);};
      intps.push_back(AdaptiveInterpolator());
      intps.at(i).initalize(func, fix.x1, fix.x2, fix.mean_error_bound); 
    }
    for (auto _ : state) {
      for (int j{0}; j<nb_intps; j++) {
        for (int i{0}; i<fix.points.size();i++) {
          intps[j].interpolate(fix.points(i));
        }
      }
    }
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"nb_points",fix.nb_points}
      });
  }



  template <class Fix>
  void BM_Hyp1f1(benchmark::State &state) {
    auto fix = Fix(state);
    double n = 10;
    double l = 10;
    double a = 0.5*(n+l+3);
    double b = l+1.5;
    auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
    for (auto _ : state) {
      for (int i{0}; i<fix.points.size();i++) {
        hyp1f1.calc(fix.points(i));
      }
    }
    state.counters.insert({
        {"nb_points",fix.nb_points}
      });
  }

  template <class Fix>
  void BM_InterpolatorHyp1f1(benchmark::State &state) {
    auto fix = Fix(state);
    auto intp{Interpolator <
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
      SearchMethod<SearchMethod_t::Hunt>
        >()};
    double n = 10;
    double l = 10;
    double a = 0.5*(n+l+3);
    double b = l+1.5;
    auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
    auto func = [&hyp1f1](double x) {return hyp1f1.calc(x);};
    intp.initalize(func, fix.x1, fix.x2, fix.mean_error_bound); 
    for (auto _ : state) {
      for (int i{0}; i<fix.points.size();i++) {
        intp.interpolate(fix.points(i));
      }
    }
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"log(mean_error)",std::log10(intp.mean_error)},
        {"log(max_error)",std::log10(intp.max_error)},
        {"nb_points",fix.nb_points},
        {"grid_size",intp.grid_rational.grid_size}
      });
  }


  template <class Fix>
  void BM_IntpSimpFun(benchmark::State &state) {
    auto fix = Fix(state);
    auto intp{Interpolator <
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
      SearchMethod<SearchMethod_t::Hunt>
        >()};
    intp.initalize(fix.func, fix.x1, fix.x2, fix.mean_error_bound); 
    for (auto _ : state) {
      for (int i{0}; i<fix.points.size();i++) {
        intp.interpolate(fix.points(i));
      }
    }
    state.counters.insert({
        {"x1",fix.x1},
        {"x2",fix.x2},
        {"log(mean_error_bound)", fix.log_mean_error_bound},
        {"log(mean_error)",std::log10(intp.mean_error)},
        {"log(max_error)",std::log10(intp.max_error)},
        //{"nb_points",fix.nb_points},
        {"grid_size",intp.grid_rational.grid_size}
      });
  }


  BENCHMARK_TEMPLATE(BM_IntpSimpFun, IntpFix<I_B>)->Apply(AllCombinationsArguments<I_B>);

  //BENCHMARK_TEMPLATE(BM_RadialContr, IntpFix<FunctionSmallDataset>)->Apply(AllCombinationsArguments<FunctionSmallDataset>);
  //BENCHMARK_TEMPLATE(BM_Hyp1f1, IntpFix<FunctionSmallDataset>)->Apply(AllCombinationsArguments<FunctionSmallDataset>);

  //BENCHMARK_TEMPLATE(BM_InterpolatorHyp1f1, IntpFix<InterpolatorSmallDataset>)->Apply(AllCombinationsArguments<InterpolatorSmallDataset>);

  //BENCHMARK_TEMPLATE(BM_InterpolatorRadialContr, IntpFix<SmallDataset>)->Apply(AllCombinationsArguments<FunctionSmallDataset>);





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
