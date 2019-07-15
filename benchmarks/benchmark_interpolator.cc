#include <benchmark/benchmark.h>
#include <functional>
#include <map>
#include <iostream>

#include "math/interpolator.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "json.hpp"



using namespace rascal::math;

namespace rascal {
  namespace internal {
    //////////
    // Issues //
    // //////
    // Currently there is not way to print out the results in scientific notation which make it harder to interpret the results. There is a pull request takling this issue https://github.com/google/benchmark/pull/821 but it is not yet merged
    
    // TODO(alex) To have a credible time benchmarks for the initialization function.
    // distance function caluclates for (radial_max, angular_max+1) points, this must be incuded
    // in the intp, maybe make extra class MatrixIntp

    //test interpoltar matrix
    //1.tests Interpoltaro for gauss function compares only initialization speed


  // Helper function to create all combinations from two arguments, because google only offers combinations of multiple Ranges with exponential incrementation. It is a modified copy of the Ranges function in google benchmark https://github.com/google/benchmark/blob/master/src/benchmark_register.cc#L300-L332
  //
  std::vector<std::vector<int64_t>> combinations(std::vector<std::vector<int64_t>> arglists) {
    //arglists = {{0,1},{0,1,2}};//(ranges.size());
   
    std::size_t total = 1;
    for (std::size_t i = 0; i < arglists.size(); i++) {
      total *= arglists[i].size();
    }
    std::vector<std::vector<int64_t>> args_;
    
    std::vector<std::size_t> ctr(arglists.size(), 0);

    for (std::size_t i = 0; i < total; i++) {
      std::vector<int64_t> tmp;
      tmp.reserve(arglists.size());

      for (std::size_t j = 0; j < arglists.size(); j++) {
        tmp.push_back(arglists[j].at(ctr[j]));
      }

      args_.push_back(std::move(tmp));

      for (std::size_t j = 0; j < arglists.size(); j++) {
        if (ctr[j] + 1 < arglists[j].size()) {
          ++ctr[j];
          break;
        }
        ctr[j] = 0;
      }
    }
    return args_;
  }

  // Because Ranges jumps in exponential
  template <class Dataset>
  void AllCombinationsArguments(benchmark::internal::Benchmark* b) {
    json data = Dataset::data();
    std::vector<std::vector<int64_t>> nbs_args;
    nbs_args.resize(data.size());
    // We have to cound backwards to keep the order in the benchmark as in the json string, because the function creating the combinations, creates them backwards
    size_t i{0};
    for (json::iterator it = data.begin(); it != data.end(); ++it) {
      std::cout << it.key() << "=" << it.value().size() << std::endl;
      std::vector<int64_t> args(it.value().size()) ; // vector with 100 ints.
      std::iota (std::begin(args), std::end(args), 0); // Fill with 0, 1, ..., 99.
      nbs_args.at(i) = args;
      i++;
    }
    auto args_comb = combinations(nbs_args);
    for (auto it = args_comb.begin(); it != args_comb.end(); ++it) {
      std::vector<int64_t> ve = *it;
      for (auto it2 = ve.begin(); it2 != ve.end(); ++it2) {
        std::cout << *it2 << " ";
      }
      std::cout << std::endl;
      b->Args(*it);
    }
    std::cout << std::endl;
  }




    static void BM_RadialContr(benchmark::State& state) {
      int max_radial{20};
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
      math::Vector_t points = math::Vector_t::LinSpaced(100,0,5);
      math::Vector_t pointsr = math::Vector_t::Zero(100);
      for (auto _ : state) {
        for (int i{0}; i<points.size();i++) {
          radial_contr.compute_contribution<AtomicSmearingType::Constant>(points(i), 0.5)(0,17);
        }
      }
    }
    BENCHMARK(BM_RadialContr);

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


 //class BaseFixture : public ::benchmark::Fixture {
 //  SetUp() {
 //    std::vector<size_t> sizes;
 //  }
 //  json data;
 //}
 //https://github.com/google/benchmark/issues/387
 class MyFixture : public ::benchmark::Fixture {
   public:
    //if needed to initialize reading data
    //void SetUp(const ::benchmark::State& state) {}
    //void TearDown(const ::benchmark::State& state) {}
    
    const std::vector<double> x1s{0};
    const std::vector<double> x2s{5};
    const std::vector<double> mean_error_bounds{1e-1,1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7,1e-8,1e-9, 1e-10};
    const std::vector<size_t> nbs_points{2000};
  };
  

  //This is a test which is executed one time if one or more BENCHMARK_DEFINE_F(MyFixture, ....) is used. For repeated BENCHMARK_DEFINE_F it still is only executed one time
  //BENCHMARK_F(MyFixture, Foo)(benchmark::State &st) {
  //  for (auto _ : st) {
  //    for (int i{0};i<50;i++){ if(i==data){data++;}}
  //  }
  //}
  //BENCHMARK_DEFINE_F(MyFixture, Intp)(benchmark::State &state) {
  //  data.get<double>("x0").at(state.range(0));
  //  double x1{x1s.at(state.range(0))};
  //  double x2{x2s.at(state.range(1))};
  //  double mean_error_bound{mean_error_bounds.at(state.range(2))};
  //  size_t nb_points{nbs_points.at(state.range(3))};

  //  return_sizes() {

  //  }
  //  
  //  {"x1","{0,1}"}
  //  auto intp{Interpolator <
  //    InterpolationMethod<InterpolationMethod_t::CubicSpline>,
  //    GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
  //    SearchMethod<SearchMethod_t::Hunt>
  //      >()};
  //  math::Vector_t points = math::Vector_t::LinSpaced(nb_points,x1,x2);
  //  auto exp_func{[](double x) {return std::exp(x);}};
  //  intp.initalize(exp_func, x1,x2, mean_error_bound); 
  //  for (auto _ : state) {
  //    for (int i{0}; i<points.size();i++) {
  //      intp.interpolate(points(i));
  //    }
  //  }
  //  state.counters.insert({{"x1",x1},{"x2",x2}, {"log_mean_error_bound", std::log10(mean_error_bound)},{"nb_points",nb_points},{"GridSize",intp.grid_rational.grid_size}});
  //}

  //BENCHMARK_DEFINE_F(MyFixture, FooTest2)(benchmark::State &st) {
  //  for (auto _ : st) {
  //    for (int i{0};i<100;i++){ if(i==data){data++;}}
  //  }
  //}
  
  // The data fixture contains the data structure for the benchmarks in a static function
  
  //class InterpolatorDataFixture : public DataFixture {
  // public:
  //  // static json has to be captured in function because it is a non-literal
  //  static json data() {
  //    return {
  //      {"x1x2", x1x2s},
  //      {"mean_error_bound", mean_error_bounds}
  //      };
  //  }

  // private: 
  //  static const std::vector<double> x1x2s;
  //  static const std::vector<double> mean_error_bounds;
  //};
  //
  //
  
  // GoogleBenchmark has the functionality to add all combinations between multiple Ranges of arguments. If I want to check my function all combinations Ra 
  
  class BaseInterpolatorDataset  {
   public:
    enum class SupportedFunc {
      Gauss,
      Hyp1f1
    };
  };

  /* To make the data available in an Argument function we create a static class with all the data. To avoid separate declaration and definitions for static member variables we use functions. This explained more in detail in https://stackoverflow.com/a/17057121/10329403 .
   */

  class Hyp1f1Dataset : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,3)}},
         {"log_mean_error_bounds", {-3}},
         {"func_names", {SupportedFunc::Gauss}},
         {"nbs_points", {1000,5000,10000}}
         };
     }
  };
  class InterpolatorSmallDataset : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,3)}},
         {"log_mean_error_bounds", {-3,-5,-8}},
         {"func_names", {SupportedFunc::Gauss}},
         {"nbs_points", {1000,5000,10000}}
         };
     }
  };
  class InterpolatorBigDataset : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,3),std::make_pair(0,5),std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-3,-5,-8,-10,-12}},
         {"func_names", {SupportedFunc::Gauss}},
         //{"func", {"gauss","double_gauss","inverted_double_gauss","hyp1f1", "orthonormalized_hyp1f1"}}
         {"nbs_points", {1000,5000,10000}}
         };
     }
  };
  
  template<class Dataset>
  class BaseFixture {
   protected:
    BaseFixture() {
      json data = Dataset::data();
      this->build_state_range_index_to_key(data);
    }
    // We have do this with an iterator of the json file because it does not
    // perserve the insert order.
    void build_state_range_index_to_key(json & data) {
      //state_range_index_to_key.resize(data.size()); // TODO(alex) change to capacity
      int64_t i{0};
      for (json::iterator it = data.begin(); it != data.end(); ++it) {
        this->state_range_index_to_key.insert(std::pair<std::string, int64_t>(it.key(),i));
        i++;
      }
    }
    
    template <typename T>
    T lookup(json & data, std::string name, const ::benchmark::State& state) {
      return data[name].at(state.range(this->state_range_index_to_key[name])).get<T>();
    }

    std::map<std::string,int64_t> state_range_index_to_key;
  };

  template<class Dataset>
  class InterpolatorFixture : public BaseFixture<Dataset> {
   public:
    using Parent = BaseFixture<Dataset>;
    using SupportedFunc = typename Dataset::SupportedFunc;
    
    InterpolatorFixture(const ::benchmark::State& state) : Parent() {
      // TODO(alex) make these two lines to a general SetUp
      json data = Dataset::data();
      //Parent::SetUp(data);
      //this->build_state_range_index_to_key(data);

      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->log_mean_error_bound = this->template lookup<int>(data, "log_mean_error_bounds", state);
      this->mean_error_bound = std::pow(10,this->log_mean_error_bound);
      auto func_name = this->template lookup<SupportedFunc>(data, "func_names", state);
      this->set_function(func_name);
      this->nb_points = this->template lookup<size_t>(data, "nbs_points", state);
      this->points = math::Vector_t::LinSpaced(this->nb_points,this->x1,this->x2);
    }

    double x1{0};
    double x2{0};
    int log_mean_error_bound{0};
    double mean_error_bound{0};
    std::function<double(double)> func{};
    size_t nb_points{0};
    math::Vector_t points{};

   private:
    void set_function(SupportedFunc name) {
      switch(name) {
        case SupportedFunc::Gauss:
          this->func = [](double x) {return std::exp(-std::pow((x-1)/0.5,2)/2);};
          break;
        //case SupportedFunc::Hyp1f1: {
        //  double n = 10;
        //  double l = 10;
        //  double a = 0.5*(n+l+3);
        //  double b = l+1.5;
        //  // check if this does not result into segmentation fault
        //  //this->hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
        //  //this->func = [this->hyp1f1](double x) {return hyp1f1.calc(x);};
        //  break;
        //}
        default:
          this->func = [](double x) {return x;};
          break;
      }
    }
  };

  BENCHMARK_DEFINE_F(MyFixture, BenchTests)(benchmark::State &state) {
    //std::cout << state.range(0) << std::endl;
    for (auto _ : state) {
    }
  }

  //BENCHMARK_DEFINE_F(InterpolatorDataFixture, Hyp1f1)(benchmark::State &state) {
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
  BENCHMARK_TEMPLATE(BM_Hyp1f1, InterpolatorFixture<Hyp1f1Dataset>)->Apply(AllCombinationsArguments<Hyp1f1Dataset>);

  //BENCHMARK_DEFINE_F(IntWithDat, BenchTests2)(benchmark::State &state) {
  //    //this->log_mean_error_bound = data["log_mean_error_bounds"].at(state.range(0)).get<std::pair<double,double>>();
  //  //double x1 = x1x2s().at(state.range(0));
  //  auto intp{Interpolator <
  //    InterpolationMethod<InterpolationMethod_t::CubicSpline>,
  //    GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
  //    SearchMethod<SearchMethod_t::Hunt>
  //      >()};
  //  state.counters.insert({
  //      {"x1",this->x1},
  //      {"x2",this->x2},
  //      {"log(mean_error_bound)", this->log_mean_error_bound},
  //      {"log(mean_error)",std::log10(intp.mean_error)},
  //      {"log(max_error)",std::log10(intp.max_error)},
  //      {"nb_points",this->nb_points},
  //      {"grid_size",intp.grid_rational.grid_size}
  //    });
  //}

  template <class Fix>
  void BM_Interpolator(benchmark::State &state) {
      //this->log_mean_error_bound = data["log_mean_error_bounds"].at(state.range(0)).get<std::pair<double,double>>();
    //double x1 = x1x2s().at(state.range(0));
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
  BENCHMARK_TEMPLATE(BM_Interpolator, InterpolatorFixture<InterpolatorSmallDataset>)->Apply(AllCombinationsArguments<InterpolatorSmallDataset>);

  ////BENCHMARK_REGISTER_F(InterpolatorDataFixture, BenchTests2)->Apply(AllCombinationsArguments<InterpolatorDataFixture>);

  //using IntWithDat = InterpolatorDataFixture<DataFixture>;
  //BENCHMARK_REGISTER_F(IntWithDat, BenchTests3)->Apply(AllCombinationsArguments<IntWithDat>);



  //// for one argument DenseRange can be used
  ////BENCHMARK_REGISTER_F(MyFixture, Intp)->Appyl(IteratorArguments<MyFixture>)
  ////  Ranges({{0,0},{0,0},{0,10},{0,0}});
  ////BENCHMARK_REGISTER_F(MyFixture, FooTest2);

  //  //parameters x1 x2 precision, nb_test_points
  //  //additional UserCounter
  //  //grid_size,
  //  
  //  // mean_error_bound
  //  // max_error
  //  void BM_Intp(benchmark::State& state, double x1, double x2, int log_mean_error_bound, int nb_points) {
  //    auto intp{Interpolator <
  //      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
  //      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
  //      SearchMethod<SearchMethod_t::Hunt>
  //        >()};
  //    math::Vector_t points = math::Vector_t::LinSpaced(nb_points,x1,x2);
  //    auto exp_func{[](double x) {return std::exp(x);}};
  //    double mean_error_bound{std::pow(10,log_mean_error_bound)};
  //    intp.initalize(exp_func, x1,x2, mean_error_bound); 
  //    for (auto _ : state) {
  //      for (int i{0}; i<points.size();i++) {
  //        intp.interpolate(points(i));
  //      }
  //    }
  //    state.counters.insert({{"x1",x1},{"x2",x2},{"GridSize",intp.grid_rational.grid_size}});
  //  }
  //  BENCHMARK_CAPTURE(BM_Intp, basic_test, 0,5,1e-3,1000);

  //  void BM_RadialContrIntp(benchmark::State& state, double x1, double x2, double precision, int nb_points) {
  //  auto intp{Interpolator <
  //    InterpolationMethod<InterpolationMethod_t::CubicSpline>,
  //    GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
  //    SearchMethod<SearchMethod_t::Hunt>
  //      >()};

  //    const json fc_hypers{
  //         {"type", "Constant"},
  //         {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
  //        };
  //    const json hypers{{"gaussian_density", fc_hypers},
  //              {"max_radial", 20},
  //              {"max_angular", 19},
  //              {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
  //    };
  //    auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
  //    std::function<double(double)> func = [&radial_contr](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x,0.5)(0,17);};
  //    math::Vector_t points = math::Vector_t::LinSpaced(nb_points,x1,x2);
  //    intp.initalize(func, x1,x2, std::pow(10,-precision)); 
  //    for (auto _ : state) {
  //      for (int i{0}; i<points.size();i++) {
  //        intp.interpolate(points(i));
  //      }
  //    }
  //    state.counters.insert({{"x1",x1},{"x2",x2},{"-log(precision)",precision},{"nb_points",nb_points},{"GridSize",intp.grid_rational.grid_size}});
  //  }
  //  //BENCHMARK_CAPTURE(BM_RadialContrIntp, basic_test, 0,5,5,1000);
  //    //->Args({0.,5.,1e-5,1000}); // does not work because Argas only accepts integers
  //    //->Ranges({{0,0},{5,5},{1e-2, 1e-12},{100,3000}});
  //  //{1e-2,1e-4,1e-6,1e-8,1e-10,1e-12}
  //  // 100 times slower than with precompution that is why we skip
  //  


  //  void BM_RadialContrIntpS(benchmark::State& state, double x1, double x2, double precision, int nb_points) {
  //  auto intp{Interpolator <
  //    InterpolationMethod<InterpolationMethod_t::CubicSpline>,
  //    GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>,
  //    SearchMethod<SearchMethod_t::Hunt>
  //      >()};
  //    std::function<double(double)> func = [](double x) {return std::exp(-std::pow((x-1)/.5,2)/2);};
  //    math::Vector_t points = math::Vector_t::LinSpaced(nb_points,x1,x2);
  //    intp.initalize(func, x1,x2, precision); 
  //    for (auto _ : state) {
  //      for (int i{0}; i<points.size();i++) {
  //        intp.interpolate(points(i));
  //      }
  //    }
  //    state.counters.insert({{"x1",x1},{"x2",x2},{"precision",precision},{"nb_points",nb_points},{"GridSize",intp.grid_rational.grid_size}});
  //  }
  //  
  //  BENCHMARK_CAPTURE(BM_RadialContrIntpS, basic_test, 0,5,1e-5,1000);
    //static double radial_contr_function_generator(int n, int l, double r) {

    //  int max_radial{20};
    //  int max_angular{max_radial-1};
    //  json fc_hypers{
    //       {"type", "Constant"},
    //       {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
    //      };
    //  json hypers{{"gaussian_density", fc_hypers},
    //            {"max_radial", max_radial},
    //            {"max_angular", max_angular},
    //            {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
    //  };
    //  auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
    //  return radial_contr.compute_contribution<AtomicSmearingType::Constant>(r, 0.5)(n,l);
    //}
    //static void BM_RadialContrFunc(benchmark::State& state) {
    //  std::function<double(double)> func{std::bind(radial_contr_function_generator, 0, 17, std::placeholders::_1)};
    //  math::Vector_t points = math::Vector_t::LinSpaced(3000,0,5);
    //  for (auto _ : state) {
    //    for (int i{0}; i<points.size();i++) {
    //      func(points(i));
    //    }
    //  }
    //}
    // Register the function as a benchmark
    //BENCHMARK(BM_RadialContrFunc);
  }
}

BENCHMARK_MAIN();
