#include <functional>
#include <map>
#include <iostream>

#include "math/interpolator.hh"
#include "benchmarks.hh"  
#include "representations/representation_manager_spherical_expansion.hh"
#include "json.hpp"

using namespace rascal::math;

namespace rascal {
  namespace internal {
  int SEED = 1597463007; //0x5f3759df

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
         {"nbs_points", {1e3,1e4,1e5,1e6}},
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
         {"log_mean_error_bounds", {-10}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_points", {1e3,1e4,1e5,1e6}},
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
      //json data = Dataset::data();
      //auto range = this->template lookup<std::pair<double,double>>(data, "ranges", 0); 
      //this->x1 = std::get<0>(range);
      //this->x2 = std::get<1>(range);
      //this->log_mean_error_bound = this->template lookup<int>(data, "log_mean_error_bounds", 0);
      //this->mean_error_bound = std::pow(10,this->log_mean_error_bound);
      //auto func_name = this->template lookup<SupportedFunc>(data, "func_names", 0);
      //this->set_function(func_name);
      //this->intp.initalize(this->func, this->x1, this->x2, this->mean_error_bound); 
      //this->init_points();
      //this->init = true;
    }

    void init_points() {
      const int new_points = 1000000;
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

    // TODO(alex) only way to prevent multiple initialization but still change of parameters when they change is to have a has_changed function which checks if parameters have changed
    // To split parameters which do not change the interpolator from nb_interations, I can include a has_interpolator_parameters_changed
    
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
        std::cout << "Init ressources " << this->log_mean_error_bound << std::endl;
      }
      this->nb_points = this->template lookup<size_t>(data, "nbs_points", state);
      this->random = this->template lookup<bool>(data, "random", state); 
      //std::cout << "Finish points" << std::endl;
    }
    void TearDown() {
      init = false;
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
