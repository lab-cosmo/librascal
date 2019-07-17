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
  int SEED = 1597463007; //0x5f3759df
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
      //std::cout << it.key() << "=" << it.value().size() << std::endl;
      std::vector<int64_t> args(it.value().size()) ; // vector with 100 ints.
      std::iota (std::begin(args), std::end(args), 0); // Fill with 0, 1, ..., 99.
      nbs_args.at(i) = args;
      i++;
    }
    auto args_comb = combinations(nbs_args);
    for (auto it = args_comb.begin(); it != args_comb.end(); ++it) {
      std::vector<int64_t> ve = *it;
      for (auto it2 = ve.begin(); it2 != ve.end(); ++it2) {
        //std::cout << *it2 << " ";
      }
      //std::cout << std::endl;
      b->Args(*it);
    }
    //std::cout << std::endl;
  }

  class BaseInterpolatorDataset  {
   public:
    // TODO(alex) better naming
    enum class SupportedFunc {
      Gaussian,
      TwoGaussians,
      SinLikeGaussian
    };
  };

  /* To make the data available in an Argument function we create a static class with all the data. To avoid separate declaration and definitions for static member variables we use functions. This explained more in detail in https://stackoverflow.com/a/17057121/10329403 .
   */

  // ClassFunctionSmallDataset, used for functions capsulated 
  class CF_S : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-3}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {1e6}},
         {"random", {true}}
         };
     }
  };
  class CF_B : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-10}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {1e6}},
         {"random", {true}}
         };
     }
  };
  class CFI_S : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-3}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {1e6}},
         {"random", {true}}
         };
     }
  };
  class CFI_B : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-9}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {1e6}},
         {"random", {true}}
         };
     }
  };

  class I_S : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-10}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {100,1000,10000,100000}},
         {"random", {false}}
         };
     }
  };
  class I_B : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,3),std::make_pair(0,5),std::make_pair(0,8)}},
         {"log_mean_error_bounds", {-3,-5,-8}},
         {"func_names", {SupportedFunc::Gaussian, SupportedFunc::TwoGaussians,SupportedFunc::SinLikeGaussian}},
         {"nbs_points", {1e8}},
         {"random", {true}}
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
      return this->template lookup<T>(data, name, state.range(this->state_range_index_to_key[name])); 
    }

    template <typename T>
    T lookup(json & data, std::string name, const int64_t index) {
      return data[name].at(index).get<T>();
    }

    std::map<std::string,int64_t> state_range_index_to_key;
  };

  template<class Dataset>
  class IntpFix : public BaseFixture<Dataset> {
   public:
    using Parent = BaseFixture<Dataset>;
    using SupportedFunc = typename Dataset::SupportedFunc;
    using Interpolator_t = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::HeapBased>,
      SearchMethod<SearchMethod_t::AStarUniform>
        >;
    
    // global init
    IntpFix<Dataset>() : Parent() {
      json data = Dataset::data();
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", 0); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->log_mean_error_bound = this->template lookup<int>(data, "log_mean_error_bounds", 0);
      this->mean_error_bound = std::pow(10,this->log_mean_error_bound);
      auto func_name = this->template lookup<SupportedFunc>(data, "func_names", 0);
      this->set_function(func_name);
      this->intp.initalize(this->func, this->x1, this->x2, this->mean_error_bound); 
      this->init_points();
      this->init =true;
    }

    void init_points(){
      const int new_points = 100000;
      if (this->random) {
        srand(SEED);
        math::Vector_t points_tmp  = math::Vector_t::LinSpaced(new_points , this->x1,this->x2);
        this->points = math::Vector_t::Zero(new_points);
        for (int i{0};i<this->points.size();i++) {
          this->points(i) = points_tmp(rand() % new_points);
        }
      } else {
        this->points = math::Vector_t::LinSpaced(new_points,this->x1,this->x2);
      }
    }

    // local init TODO(alex) remove
    IntpFix<Dataset>(const ::benchmark::State& state) : Parent() {
      this->SetUp(state);
    }

    // local init
    void SetUp(const ::benchmark::State& state) {
      json data = Dataset::data();
      if (not(this->init)) {
        // TODO(alex) make these two lines to a general SetUp
        //Parent::SetUp(data);
        //this->build_state_range_index_to_key(data);

        auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
        this->x1 = std::get<0>(range);
        this->x2 = std::get<1>(range);
        this->log_mean_error_bound = this->template lookup<int>(data, "log_mean_error_bounds", state);
        this->mean_error_bound = std::pow(10,this->log_mean_error_bound);
        auto func_name = this->template lookup<SupportedFunc>(data, "func_names", state);
        this->set_function(func_name);
        this->intp.initalize(this->func, this->x1, this->x2, this->mean_error_bound); 
        this->init_points();
        this->init =true;
      }
      std::cout << "Init ressources" << this->log_mean_error_bound << std::endl;
      this->nb_points = this->template lookup<size_t>(data, "nbs_points", state);
      this->random = this->template lookup<bool>(data, "random", state); 
      //std::cout << "Finish points" << std::endl;
    }
    void TearDown() {
      this->init = false;
    }

    bool init;
    Interpolator_t intp;
    double x1{0};
    double x2{0};
    int log_mean_error_bound{0};
    double mean_error_bound{0};
    std::function<double(double)> func{};
    size_t nb_points{0};
    bool random;
    math::Vector_t points{};

   private:
    void set_function(SupportedFunc name) {
      switch(name) {
        case SupportedFunc::Gaussian:
          this->func = [](double x) {return std::exp(-std::pow((x-1)/0.5,2)/2);};
          break;
        case SupportedFunc::TwoGaussians:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) + std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        case SupportedFunc::SinLikeGaussian:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) - std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        default:
          this->func = [](double x) {return x;};
          break;
      }
    }
  };

  }
}
