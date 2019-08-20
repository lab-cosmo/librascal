#include <functional>
#include <map>
#include <iostream>

#include "math/interpolator.hh"
#include "benchmarks.hh"  
#include "representations/representation_manager_spherical_expansion.hh"
#include "json.hpp"

using namespace rascal::math;

// TODO(alex) naming to BFixture to prevent collision

namespace rascal {
  namespace internal {
  int SEED = 1597463007; //0x5f3759df

  class BaseInterpolatorDataset  {
   public:
    // TODO(alex) better naming
    enum class SupportedFunc {
      Identity,
      Gaussian,
      TwoGaussians, // TODO(alex) better naming
      SinLikeGaussian, // TODO(alex) better naming
      Hyp1f1, 
    };
    enum class SupportedVecFunc {
      RadialContribution
    };
  };

  /* To make the data available in an Argument function we create a static class with all the data. To avoid separate declaration and definitions for static member variables we use functions. This explained more in detail in https://stackoverflow.com/a/17057121/10329403 .
   */


  class RadConDataset : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         //{"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-10}},
         {"func_names", {SupportedVecFunc::RadialContribution}},
         {"max_radial", {3}},
         {"random", {true}}
         };
     }
  };

  class RadConDataset1 : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         //{"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-10}},
         {"func_names", {SupportedVecFunc::RadialContribution}},
         {"max_radial", {3}},
         {"random", {true}}
         };
     }
  };

  class RadConDataset2 : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         //{"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-10}},
         {"func_names", {SupportedVecFunc::RadialContribution}},
         {"max_radial", {5}},
         {"random", {true}}
         };
     }
  };
  class RadConDataset3 : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         //{"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-10}},
         {"func_names", {SupportedVecFunc::RadialContribution}},
         {"max_radial", {8}},
         {"random", {true}}
         };
     }
  };



  class Hyp1f1Dataset : public BaseInterpolatorDataset {
    public:
     using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
     static const json data() {
       return {
         {"nbs_iterations", {1e3,1e4,1e5,1e6}},
         {"ranges", {std::make_pair(0,16)}},
         {"log_error_bounds", {-8}},
         {"func_names", {SupportedFunc::Hyp1f1}},
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

    // Could be moved to base class with virtual classes and shared by vectorized and scalar
    void SetUp(const ::benchmark::State& state) {
      const json data = Dataset::data();
      // Because in the two initialization processes share parameters of the json string, therefore we check the change of parameters before anything is initialized
      bool interpolator_parameters_changed{this->have_interpolator_parameters_changed(state, data)};
      bool ref_points_parameters_changed{this->have_ref_points_parameters_changed(state, data)};

      if (not(this->initialized) || interpolator_parameters_changed) {
        this->init_interpolator(state, data);
      }
      if (not(this->initialized) || ref_points_parameters_changed) {
        this->init_ref_points(state, data);      
      }
      this->nb_iterations = this->template lookup<size_t>(data, "nbs_iterations", state);
      this->initialized = true;
    }
    
    bool initialized{false};
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

   protected:
    // could be moved to a base class
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
        math::Vector_t points_tmp  = math::Vector_t::LinSpaced(this->nb_ref_points, this->x1,this->x2);
        this->ref_points = math::Vector_t::Zero(this->nb_ref_points);
        for (int i{0};i<this->ref_points.size();i++) {
          this->ref_points(i) = points_tmp(rand() % this->nb_ref_points);
        }
      } else {
        this->ref_points = math::Vector_t::LinSpaced(this->nb_ref_points,this->x1,this->x2);
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
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->log_error_bound = this->template lookup<int>(data, "log_error_bounds", state);
      this->func_name = this->template lookup<SupportedFunc>(data, "func_names", state);
        
      this->error_bound = std::pow(10,this->log_error_bound);
      this->init_function();
      this->intp.initialize(this->func, this->x1, this->x2, this->error_bound); 
    }

    void init_hyp1f1_function() {
      double n = 10;
      double l = 9;
      double a = 0.5*(n+l+3);
      double b = l+1.5;
      auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
      this->func = [=](double x) mutable {return hyp1f1.calc(x);};
    }

    void init_function() {
      switch(this->func_name) {
        case SupportedFunc::Identity:
          this->func = [](double x) {return x;};
          break;
        case SupportedFunc::Gaussian:
          this->func = [](double x) {return std::exp(-std::pow((x-1)/0.5,2)/2);};
          break;
        // TODO(alex) remove these two
        case SupportedFunc::TwoGaussians:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) + std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        case SupportedFunc::SinLikeGaussian:
          this->func = [](double x) {return (std::exp(-std::pow((x-1)/0.5,2)/2) - std::exp(-std::pow((x-3)/0.5,2)/2))/2;};
          break;
        case SupportedFunc::Hyp1f1:
          this->init_hyp1f1_function();
          break;
        default:
          this->func = [](double x) {return x;};
          break;
      }
    }
  };

  template<class Dataset>
  class InterpolatorVectorizedFixture : public InterpolatorFixture<Dataset> {
   public:
    using Parent = InterpolatorFixture<Dataset>;
    using SupportedVecFunc = typename Dataset::SupportedVecFunc;
    using Interpolator_t = InterpolatorVectorized<
      InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
      SearchMethod<SearchMethod_t::Uniform>
        >;

    InterpolatorVectorizedFixture<Dataset>() : Parent() {}

    void SetUp(const ::benchmark::State& state) {
      const json data = Dataset::data();
      // Because in the two initialization processes share parameters of the json string, therefore we check the change of parameters before anything is initialized
      bool interpolator_parameters_changed{this->have_interpolator_parameters_changed(state, data)};
      bool ref_points_parameters_changed{this->have_ref_points_parameters_changed(state, data)};

      if (not(this->initialized) || interpolator_parameters_changed) {
        this->init_interpolator(state, data);
      }
      if (not(this->initialized) || ref_points_parameters_changed) {
        this->init_ref_points(state, data);      
      }
      this->nb_iterations = this->template lookup<size_t>(data, "nbs_iterations", state);
      this->initialized = true;
    }

    Interpolator_t intp;  
    std::function<math::Matrix_t(double)> func;
    SupportedVecFunc func_name;
    int max_radial{0};
    RadialContribution<RadialBasisType::GTO> radial_contr;

   protected:
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
        math::Vector_t points_tmp  = math::Vector_t::LinSpaced(this->nb_ref_points, this->x1,this->x2);
        this->ref_points = math::Vector_t::Zero(this->nb_ref_points);
        for (int i{0};i<this->ref_points.size();i++) {
          this->ref_points(i) = points_tmp(rand() % this->nb_ref_points);
        }
      } else {
        this->ref_points = math::Vector_t::LinSpaced(this->nb_ref_points,this->x1,this->x2);
      }
    }

    bool have_interpolator_parameters_changed(const ::benchmark::State& state, const json & data) const {
      bool have_scalar_parameters_changed{Parent::have_interpolator_parameters_changed(state, data)};    
      int new_max_radial = this->template lookup<int>(data, "max_radial", state);
      return (have_scalar_parameters_changed || new_max_radial != this->max_radial); 

    }

    void init_interpolator(const ::benchmark::State& state, const json & data) {
      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->log_error_bound = this->template lookup<int>(data, "log_error_bounds", state);
      this->func_name = this->template lookup<SupportedVecFunc>(data, "func_names", state);
        
      this->error_bound = std::pow(10,this->log_error_bound);
      this->init_function(state, data);
      this->intp.initialize(this->func, this->x1, this->x2, this->error_bound, 10000000, 5, true) ; 
    }

    void init_radial_contribution_function(const ::benchmark::State& state, const json & data) {
      this->max_radial = this->template lookup<int>(data, "max_radial", state);
      int max_angular{this->max_radial};
      json fc_hypers{
           {"type", "Constant"},
           {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
          };
      json hypers{{"gaussian_density", fc_hypers},
                {"max_radial", max_radial},
                {"max_angular", max_angular},
                {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
        };
      // we cannot copy radial contribution to the lambda function, because the copying has been disabled
      this->radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
      this->func = [&](double x) mutable {return this->radial_contr.compute_contribution<AtomicSmearingType::Constant>(x,0.5);};
    }

    void init_function(const ::benchmark::State& state, const json & data) {
      switch(this->func_name) {
        // This case uses the RadialContribution class as comparisment and can therefore directly be computed for different distances
        case SupportedVecFunc::RadialContribution:
         this->init_radial_contribution_function(state, data);
         break;
      }
    }
  };

//  // This class uses the RepresentationManager class as basis to call the RadialContribution. Therefore benchmarks for different atomic structures can be done.
//  template<class Dataset>
//  class RepresentationManagerBFixture : public InterpolatorFixture<Dataset> {
//
//    using Representation_t = RepresentationManagerSphericalExpansion<
//        AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;
//    using Manager_t = 
//
//    InterpolatorRepresentationManagerFixture<Dataset>() : Parent() {}
//
//    Manager_t manager{};
//    Representation_t representation{};
//    void init_repr_calc_function(const ::benchmark::State& state, const json & data) {
//      this->max_radial = this->template lookup<int>(data, "max_radial", state);
//      int max_angular{this->max_radial};
//      // TODO(alex) put these two into the data class
//      std::string filename{"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
//      double cutoff{8.};      
//      // make structure manager
//      json hypers{{"max_radial", this->max_radial},
//                  {"max_angular", max_angular},
//                  {"soap_type", "PowerSpectrum"},
//                  {"normalize", true},
//                  {"compute_gradients", true}};
//
//      json fc_hypers{{"type", "Cosine"},
//                     {"cutoff", {{"value", cutoff}, {"unit", "AA"}}},
//                     {"smooth_width", {{"value", 0.}, {"unit", "AA"}}}};
//      json sigma_hypers{{"type", "Constant"},
//                        {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};
//
//      hypers["cutoff_function"] = fc_hypers;
//      hypers["gaussian_density"] = sigma_hypers;
//
//      hypers["radial_contribution"] = {{"type", "GTO"}};
//
//      json structure{};
//      json adaptors;
//      json ad1{{"name", "AdaptorNeighbourList"},
//               {"initialization_arguments",
//                {{"cutoff", cutoff}, {"consider_ghost_neighbours", false}}}};
//      json ad2{{"name", "AdaptorStrict"},
//               {"initialization_arguments", {{"cutoff", cutoff}}}};
//      adaptors.emplace_back(ad1);
//      adaptors.emplace_back(ad2);
//
//      AtomicStructure<3> atomic_structure{};
//      atomic_structure.set_structure(filename);
//      this->manager =
//          make_structure_manager_stack<StructureManagerCenters,
//                                       AdaptorNeighbourList, AdaptorStrict>(
//          structure, adaptors);
//      // make representation manager
//      this->representation = Representation_t(manager, hypers);
//    }
//  };
//  
//
//  // abstract class for all fixtures using the interpolator
//  template<class Dataset>
//  class InterpolatorInterfaceBF : public BaseFixture<Dataset> {
//   public:
//    using Parent = BaseFixture<Dataset>;
//    using SupportedFunc = typename Dataset::SupportedFunc;
//    using Interpolator_t = Interpolator<
//      InterpolationMethod<InterpolationMethod_t::CubicSpline>,
//      GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
//      SearchMethod<SearchMethod_t::Uniform>
//        >;
//    
//    // Could be moved to base class with virtual classes and shared by vectorized and scalar
//    void SetUp(const ::benchmark::State& state) {
//      const json data = Dataset::data();
//      // Because in the two initialization processes share parameters of the json string, therefore we check the change of parameters before anything is initialized
//      bool interpolator_parameters_changed{this->have_interpolator_parameters_changed(state, data)};
//      bool ref_points_parameters_changed{this->have_ref_points_parameters_changed(state, data)};
//
//      if (not(this->initialized) || interpolator_parameters_changed) {
//        this->init_interpolator(state, data);
//      }
//      if (not(this->initialized) || ref_points_parameters_changed) {
//        this->init_ref_points(state, data);      
//      }
//      this->nb_iterations = this->template lookup<size_t>(data, "nbs_iterations", state);
//      this->initialized = true;
//    }
//    
//    bool initialized{false};
//    Interpolator_t intp;  
//    double x1{0};
//    double x2{0};
//    int log_error_bound{0};
//    double error_bound{0};
//    size_t nb_iterations{0};
//    bool random{true};
//    const int nb_ref_points = 100000;
//    math::Vector_t ref_points{Vector_t::Zero(nb_ref_points)};
//
//
//    //SupportedFunc func_name{SupportedFunc::Identity};
//    //std::function<double(double)> func{};
//
//   protected:
//    // could be moved to a base class
//    bool have_ref_points_parameters_changed(const ::benchmark::State& state, const json & data) const {
//      bool new_random = this->template lookup<bool>(data, "random", state); 
//      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
//      double new_x1 = std::get<0>(range);
//      double new_x2 = std::get<1>(range);
//      return (new_random != this->random || new_x1 != this->x1 || new_x2 != this->x2);
//    }
//
//    // initialize the ref points 
//    void init_ref_points(const ::benchmark::State& state, const json & data) { 
//      auto range = this->template lookup<std::pair<double,double>>(data, "ranges", state); 
//      this->x1 = std::get<0>(range);
//      this->x2 = std::get<1>(range);
//      this->random = this->template lookup<bool>(data, "random", state); 
//      if (this->random) {
//        srand(SEED);
//        math::Vector_t points_tmp  = math::Vector_t::LinSpaced(this->nb_ref_points, this->x1,this->x2);
//        this->ref_points = math::Vector_t::Zero(this->nb_ref_points);
//        for (int i{0};i<this->ref_points.size();i++) {
//          this->ref_points(i) = points_tmp(rand() % this->nb_ref_points);
//        }
//      } else {
//        this->ref_points = math::Vector_t::LinSpaced(this->nb_ref_points,this->x1,this->x2);
//      }
//    }
//    virtual bool have_interpolator_parameters_changed(const ::benchmark::State& state, const json & data) const = 0;
//    virtual void init_interpolator(const ::benchmark::State& state, const json & data) = 0;
//    };
  }
}
