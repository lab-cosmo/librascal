/**
 * file  benchmark_interpolator.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief contains the data and fixtures of the interpolator benchmarks
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

#include <functional>
#include <map>
#include <iostream>

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/make_structure_manager.hh"

#include "math/interpolator.hh"
#include "benchmarks.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "json.hpp"

using namespace rascal::math;

// TODO(all) naming to BFixture to prevent collision with tests ok?

namespace rascal {
  namespace internal {
    int SEED = 1597463007;  // 0x5f3759df

    class BaseInterpolatorDataset {
     public:
      enum class SupportedFunc { Identity, Gaussian, Hyp1f1 };
      enum class SupportedVecFunc { RadialContribution };
    };

    // TODO(alex) move this to tutorial
    // Static const functions because json type cannot be static const. See
    // https://stackoverflow.com/a/17057121/10329403

    class SphericalExpansionDataset : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {
            //{"nbs_iterations", {1e3,1e4,1e5,1e6}},
            {"nbs_iterations", {1e3}},
            {"ranges", {std::make_pair(0, 16)}},  // dummy
            {"log_error_bounds", {-8}},
            {"func_names", {SupportedVecFunc::RadialContribution}},  // dummy
            {"radial_angular",
             {std::make_pair(3, 3), std::make_pair(6, 6),
              std::make_pair(9, 9)}},
            {"random", {true}},  // dummy
            {"filenames",
             {"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"}},
            {"cutoffs", {8}}};
      }
    };

    class RadConDataset : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {//{"nbs_iterations", {1e3,1e4,1e5,1e6}},
                {"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                {"ranges", {std::make_pair(0, 16)}},
                {"log_error_bounds", {-10}},
                {"func_names", {SupportedVecFunc::RadialContribution}},
                {"max_radial", {3}},
                {"random", {true}}};
      }
    };

    class RadConDataset1 : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {//{"nbs_iterations", {1e3,1e4,1e5,1e6}},
                {"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                {"ranges", {std::make_pair(0, 16)}},
                {"log_error_bounds", {-10}},
                {"func_names", {SupportedVecFunc::RadialContribution}},
                {"max_radial", {3}},
                {"random", {true}}};
      }
    };

    class RadConDataset2 : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {//{"nbs_iterations", {1e3,1e4,1e5,1e6}},
                {"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                {"ranges", {std::make_pair(0, 16)}},
                {"log_error_bounds", {-10}},
                {"func_names", {SupportedVecFunc::RadialContribution}},
                {"max_radial", {5}},
                {"random", {true}}};
      }
    };
    class RadConDataset3 : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {//{"nbs_iterations", {1e3,1e4,1e5,1e6}},
                {"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                {"ranges", {std::make_pair(0, 16)}},
                {"log_error_bounds", {-10}},
                {"func_names", {SupportedVecFunc::RadialContribution}},
                {"max_radial", {8}},
                {"random", {true}}};
      }
    };

    class Hyp1f1Dataset : public BaseInterpolatorDataset {
     public:
      using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
      static const json data() {
        return {{"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                {"ranges", {std::make_pair(0, 16)}},
                {"log_error_bounds", {-8}},
                {"func_names", {SupportedFunc::Hyp1f1}},
                {"random", {true}}};
      }
    };

    // TODO(alex) make abstract class which follows general pattern in tutorial

    // abstract class for all fixtures using the interpolator
    template <class Dataset>
    class InterpolatorBFixture : public BaseFixture<Dataset> {
     public:
      using Parent = BaseFixture<Dataset>;
      using SupportedFunc = typename Dataset::SupportedFunc;

      // Could be moved to base class with virtual classes and shared by
      // vectorized and scalar
      void SetUp(const ::benchmark::State & state) {
        const json data = Dataset::data();
        // Because in the two initialization processes share parameters of the
        // json string, therefore we check the change of parameters before
        // anything is initialized
        bool interpolator_parameters_changed{
            this->have_interpolator_parameters_changed(state, data)};
        bool ref_points_parameters_changed{
            this->have_ref_points_parameters_changed(state, data)};

        if (not(this->initialized) || interpolator_parameters_changed) {
          this->init_interpolator(state, data);
        }
        if (not(this->initialized) || ref_points_parameters_changed) {
          this->init_ref_points(state, data);
        }
        this->nb_iterations =
            this->template lookup<size_t>(data, "nbs_iterations", state);
        this->initialized = true;
      }

      bool initialized{false};
      double x1{0};
      double x2{0};
      int log_error_bound{0};
      double error_bound{0};
      size_t nb_iterations{0};
      bool random{true};
      const int nb_ref_points = 100000;
      math::Vector_t ref_points{Vector_t::Zero(nb_ref_points)};

      // SupportedFunc func_name{SupportedFunc::Identity};
      // std::function<double(double)> func{};

     protected:
      // could be moved to a base class
      bool have_ref_points_parameters_changed(const ::benchmark::State & state,
                                              const json & data) const {
        bool new_random = this->template lookup<bool>(data, "random", state);
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        double new_x1 = std::get<0>(range);
        double new_x2 = std::get<1>(range);
        return (new_random != this->random || new_x1 != this->x1 ||
                new_x2 != this->x2);
      }

      // initialize the ref points
      void init_ref_points(const ::benchmark::State & state,
                           const json & data) {
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        this->x1 = std::get<0>(range);
        this->x2 = std::get<1>(range);
        this->random = this->template lookup<bool>(data, "random", state);
        if (this->random) {
          srand(SEED);
          math::Vector_t points_tmp = math::Vector_t::LinSpaced(
              this->nb_ref_points, this->x1, this->x2);
          this->ref_points = math::Vector_t::Zero(this->nb_ref_points);
          for (int i{0}; i < this->ref_points.size(); i++) {
            this->ref_points(i) = points_tmp(rand() % this->nb_ref_points);
          }
        } else {
          this->ref_points = math::Vector_t::LinSpaced(this->nb_ref_points,
                                                       this->x1, this->x2);
        }
      }
      virtual bool
      have_interpolator_parameters_changed(const ::benchmark::State & state,
                                           const json & data) const = 0;
      virtual void init_interpolator(const ::benchmark::State & state,
                                     const json & data) = 0;
    };

    template <class Dataset>
    class InterpolatorScalarFixture : public InterpolatorBFixture<Dataset> {
     public:
      using Parent = InterpolatorBFixture<Dataset>;
      using SupportedFunc = typename Dataset::SupportedFunc;
      using Interpolator_t = Interpolator<
          InterpolationMethod<InterpolationMethod_t::CubicSpline>,
          GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
          SearchMethod<SearchMethod_t::Uniform>>;

      InterpolatorScalarFixture<Dataset>() : Parent() {}

      Interpolator_t intp;
      SupportedFunc func_name{SupportedFunc::Identity};
      std::function<double(double)> func{};

     protected:
      bool
      have_interpolator_parameters_changed(const ::benchmark::State & state,
                                           const json & data) const override {
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        double new_x1 = std::get<0>(range);
        double new_x2 = std::get<1>(range);
        int new_log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        auto new_func_name =
            this->template lookup<SupportedFunc>(data, "func_names", state);
        return (new_x1 != this->x1 || new_x2 != this->x2 ||
                new_log_error_bound != this->log_error_bound ||
                new_func_name != this->func_name);
      }

      void init_interpolator(const ::benchmark::State & state,
                             const json & data) override {
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        this->x1 = std::get<0>(range);
        this->x2 = std::get<1>(range);
        this->log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        this->func_name =
            this->template lookup<SupportedFunc>(data, "func_names", state);

        this->error_bound = std::pow(10, this->log_error_bound);
        this->init_function();
        this->intp.initialize(this->func, this->x1, this->x2,
                              this->error_bound);
      }

      void init_hyp1f1_function() {
        double n = 10;
        double l = 9;
        double a = 0.5 * (n + l + 3);
        double b = l + 1.5;
        auto hyp1f1 = math::Hyp1f1(a, b, 200, 1e-15);
        this->func = [=](double x) mutable { return hyp1f1.calc(x); };
      }

      void init_function() {
        switch (this->func_name) {
        case SupportedFunc::Identity:
          this->func = [](double x) { return x; };
          break;
        case SupportedFunc::Gaussian:
          this->func = [](double x) {
            return std::exp(-std::pow((x - 1) / 0.5, 2) / 2);
          };
          break;
        case SupportedFunc::Hyp1f1:
          this->init_hyp1f1_function();
          break;
        default:
          this->func = [](double x) { return x; };
          break;
        }
      }
    };

    template <class Dataset>
    class InterpolatorVectorizedFixture : public InterpolatorBFixture<Dataset> {
     public:
      using Parent = InterpolatorBFixture<Dataset>;
      using SupportedVecFunc = typename Dataset::SupportedVecFunc;
      using Interpolator_t = InterpolatorVectorized<
          InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
          GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
          SearchMethod<SearchMethod_t::Uniform>>;

      Interpolator_t intp;
      std::function<math::Matrix_t(double)> func;
      SupportedVecFunc func_name;
      int max_radial{0};
      RadialContribution<RadialBasisType::GTO> radial_contr;

     protected:
      bool have_scalar_interpolator_parameters_changed(
          const ::benchmark::State & state, const json & data) const {
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        double new_x1 = std::get<0>(range);
        double new_x2 = std::get<1>(range);
        int new_log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        auto new_func_name =
            this->template lookup<SupportedVecFunc>(data, "func_names", state);
        return (new_x1 != this->x1 || new_x2 != this->x2 ||
                new_log_error_bound != this->log_error_bound ||
                new_func_name != this->func_name);
      }

      bool
      have_interpolator_parameters_changed(const ::benchmark::State & state,
                                           const json & data) const override {
        bool have_scalar_parameters_changed{
            this->have_scalar_interpolator_parameters_changed(state, data)};

        int new_max_radial =
            this->template lookup<int>(data, "max_radial", state);
        return (have_scalar_parameters_changed ||
                new_max_radial != this->max_radial);
      }

      void init_interpolator(const ::benchmark::State & state,
                             const json & data) override {
        auto range = this->template lookup<std::pair<double, double>>(
            data, "ranges", state);
        this->x1 = std::get<0>(range);
        this->x2 = std::get<1>(range);
        this->log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        this->func_name =
            this->template lookup<SupportedVecFunc>(data, "func_names", state);

        this->error_bound = std::pow(10, this->log_error_bound);
        this->init_function(state, data);
        this->intp.initialize(this->func, this->x1, this->x2, this->error_bound,
                              10000000, 5, true);
      }

      void init_radial_contribution_function(const ::benchmark::State & state,
                                             const json & data) {
        this->max_radial =
            this->template lookup<int>(data, "max_radial", state);
        int max_angular{this->max_radial};
        json fc_hypers{{"type", "Constant"},
                       {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}};
        json hypers{
            {"gaussian_density", fc_hypers},
            {"max_radial", max_radial},
            {"max_angular", max_angular},
            {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};
        // we cannot copy radial contribution to the lambda function, because
        // the copying has been disabled
        this->radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
        this->func = [&](double x) mutable {
          return this->radial_contr
              .compute_contribution<AtomicSmearingType::Constant>(x, 0.5);
        };
      }

      void init_function(const ::benchmark::State & state, const json & data) {
        switch (this->func_name) {
        // This case uses the RadialContribution class as comparisment and can
        // therefore directly be computed for different distances
        case SupportedVecFunc::RadialContribution:
          this->init_radial_contribution_function(state, data);
          break;
        }
      }
    };

    // This class uses the RepresentationManager class as basis to call the
    // RadialContribution. Therefore benchmarks for different atomic structures
    // can be done.
    template <class Dataset>
    class SphericalExpansionBFixture : public InterpolatorBFixture<Dataset> {
     public:
      using Representation_t = RepresentationManagerSphericalExpansion<
          AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;
      using Manager_t =
          AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      using Parent = InterpolatorBFixture<Dataset>;

      SphericalExpansionBFixture<Dataset>(const bool use_interpolator,
                                          const bool compute_gradients)
          : Parent(), use_interpolator{use_interpolator},
            compute_gradients{compute_gradients} {}

      const bool use_interpolator;
      const bool compute_gradients;
      int max_radial{0};
      int max_angular{0};
      std::string filename;
      double cutoff{0};
      int nb_neighbours{0};
      ManagerPtr_t manager{};
      // to postpone initialization of representation from the time the object
      // is created, we create a list
      std::shared_ptr<Representation_t> representation_ptr;
      json hypers;

      void init_interpolator(const ::benchmark::State & state,
                             const json & data) override {
        // set parameter
        this->log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        this->error_bound = std::pow(10, this->log_error_bound);
        auto radial_angular = this->template lookup<std::pair<int, int>>(
            data, "radial_angular", state);
        this->max_radial = std::get<0>(radial_angular);
        this->max_angular = std::get<1>(radial_angular);
        this->filename =
            this->template lookup<std::string>(data, "filenames", state);
        this->cutoff = this->template lookup<double>(data, "cutoffs", state);

        // make structure manager
        json structure{};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                 {"initialization_arguments",
                  {{"cutoff", cutoff}, {"consider_ghost_neighbours", false}}}};
        json ad2{{"name", "AdaptorStrict"},
                 {"initialization_arguments", {{"cutoff", this->cutoff}}}};
        adaptors.emplace_back(ad1);
        adaptors.emplace_back(ad2);

        AtomicStructure<3> atomic_structure{};
        atomic_structure.set_structure(this->filename);
        this->manager =
            make_structure_manager_stack<StructureManagerCenters,
                                         AdaptorNeighbourList, AdaptorStrict>(
                structure, adaptors);
        this->manager->update(atomic_structure);

        this->nb_neighbours = 0;
        for (auto center : this->manager) {
          for (auto neigh : center) {
            this->nb_neighbours++;
          }
        }

        // make representation manager
        json hypers{{"max_radial", this->max_radial},
                    {"max_angular", this->max_angular},
                    {"soap_type", "PowerSpectrum"},
                    {"normalize", true},
                    {"compute_gradients", this->compute_gradients}};

        json fc_hypers{{"type", "Cosine"},
                       {"cutoff", {{"value", this->cutoff}, {"unit", "AA"}}},
                       {"smooth_width", {{"value", 0.}, {"unit", "AA"}}}};
        json sigma_hypers{{"type", "Constant"},
                          {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};

        hypers["cutoff_function"] = fc_hypers;
        hypers["gaussian_density"] = sigma_hypers;

        if (this->use_interpolator) {
          hypers["radial_contribution"] = {
              {"type", "GTO"},
              {"optimization",
               {{"type", "Spline"},
                {"accuracy", this->error_bound},
                {"range", {{"begin", 0.}, {"end", cutoff}}}}}};
        } else {
          hypers["radial_contribution"] = {{"type", "GTO"}};
        }
        this->hypers = hypers;
        this->representation_ptr =
            std::make_shared<Representation_t>(manager, hypers);
      }

      bool
      have_interpolator_parameters_changed(const ::benchmark::State & state,
                                           const json & data) const override {
        int new_log_error_bound =
            this->template lookup<int>(data, "log_error_bounds", state);
        auto radial_angular = this->template lookup<std::pair<int, int>>(
            data, "radial_angular", state);
        int new_max_radial = std::get<0>(radial_angular);
        int new_max_angular = std::get<1>(radial_angular);
        std::string new_filename =
            this->template lookup<std::string>(data, "filenames", state);
        double new_cutoff =
            this->template lookup<double>(data, "cutoffs", state);
        return (new_log_error_bound != this->log_error_bound ||
                new_max_radial != this->max_radial ||
                new_max_angular != this->max_angular ||
                new_filename != this->filename || new_cutoff != this->cutoff);
      }
    };
  }  // namespace internal
}  // namespace rascal
