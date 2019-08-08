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
      Identity,
      Gaussian,
      TwoGaussians,
      SinLikeGaussian,
      RadialContribution
    };
  };

  /* To make the data available in an Argument function we create a static class with all the data. To avoid separate declaration and definitions for static member variables we use functions. This explained more in detail in https://stackoverflow.com/a/17057121/10329403 .
   */


  class RadialContributionInterpolatorVectorized : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-8}},
         {"func_names", {SupportedFunc::RadialContribution}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"max_angular", {3}},
         {"random", {true}}
         };
     }
  };

  class CFI_B : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-8}},
         {"func_names", {SupportedFunc::Gaussian}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"max_angular", {3}},
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
         {"log_error_bounds", {-10}},
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
         {"log_error_bounds", {-3,-5,-8}},
         {"func_names", {SupportedFunc::Gaussian, SupportedFunc::TwoGaussians,SupportedFunc::SinLikeGaussian}},
         {"nbs_points", {1e8}},
         {"random", {true}}
         };
     }
  };


  template<class Dataset>
  class InterpolatorFixture : public BaseFixture<Dataset> {
   public:
    using Parent = BaseFixture<Dataset>;
    using SupportedFunc = typename Dataset::SupportedFunc;
    using Interpolator_t = Interpolator<
      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;
    
    InterpolatorFixture <Dataset>() : Parent() {}

    void SetUp(const ::benchmark::State& state) {
      const json data = Dataset::data();
      // Because in the two initialization processes share parameters of the json string, therefore we check the change of parameters before anything is initialized
      bool init_interpolator{this->have_interpolator_parameters_changed(state, data)};
      bool init_ref_points{this->have_ref_points_parameters_changed(state, data)};
      if (init_interpolator) {
        this->init_interpolator(state, data);
      }
      if (init_ref_points) {
        this->init_ref_points(state, data);
      }
      this->nb_iterations = this->template lookup<size_t>(data, "nbs_iterations", state);
    }

    Interpolator_t intp;
    double x1{0};
    double x2{0};
    int log_error_bound{0};
    double error_bound{0};
    SupportedFunc func_name{SupportedFunc::Identity};
    std::function<double(double)> func{};
    size_t nb_iterations{0};
    bool random{true};
    const int nb_ref_points = 100000;
    math::Vector_t ref_points{Vector_t::Zero(nb_ref_points)};
    //RadialContribution<RadialBasisType::GTO> 

   private:
    bool have_ref_points_parameters_changed(const ::benchmark::State& state, const json & data) const {
      bool new_random = this->template lookup<bool>(data, "random", state); 
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      double new_x1 = std::get<0>(range);
      double new_x2 = std::get<1>(range);
      return (new_random != this->random || new_x1 != this->x1 || new_x2 != this->x2);
    }

    // initialize the ref points 
    void init_ref_points(const ::benchmark::State& state, const json & data) { 
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->random = this->template lookup<bool>(data, "random", state); 
      if (this->random) {
        srand(SEED);
        math::Vector_t points_tmp  = math::Vector_t::LinSpaced(nb_ref_points , this->x1,this->x2);
        this->ref_points = math::Vector_t::Zero(nb_ref_points);
        for (int i{0};i<this->ref_points.size();i++) {
          this->ref_points(i) = points_tmp(rand() % nb_ref_points);
        }
      } else {
        this->ref_points = math::Vector_t::LinSpaced(nb_ref_points,this->x1,this->x2);
      }
    }

    bool have_interpolator_parameters_changed(const ::benchmark::State& state, const json & data) const {
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      double new_x1 = std::get<0>(range);
      double new_x2 = std::get<1>(range);
      int new_log_error_bound = this->template lookup<int>(data, "log_error_bounds", state);
      auto new_func_name = this->template lookup<SupportedFunc>(data, "func_names", state);
      return (new_x1 != this->x1 || new_x2 != this->x2 || new_log_error_bound != this->log_error_bound || new_func_name != this->func_name);      
    }

    void init_interpolator(const ::benchmark::State& state, const json & data) {
      std::cout << "Initialize interpolator" << std::endl;
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->log_error_bound = this->template lookup<int>(data, "log_error_bounds", state);
      this->func_name = this->template lookup<SupportedFunc>(data, "func_names", state);
        
      this->error_bound = std::pow(10,this->log_error_bound);
      //this->init_function(func_name);
      //this->intp.initialize(this->func, this->x1, this->x2, this->error_bound); 
    }

    void init_function(SupportedFunc name) {
      switch(name) {
        case SupportedFunc::Identity:
          this->func = [](double x) {return x;};
          break;
        case SupportedFunc::Gaussian:
          this->func = [](double x) {return std::exp(-std::pow((x-1)/0.5,2)/2);};
          break;
        case SupportedFunc::TwoGaussians:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) + std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        case SupportedFunc::SinLikeGaussian:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) - std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        case SupportedFunc::RadialContribution:
          //init_radial_contribution(param);
            
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
