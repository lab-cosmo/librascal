/**
 * @file  performance/benchmarks/benchmark_interpolator.hh
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

#ifndef PERFORMANCE_BENCHMARKS_BENCHMARK_INTERPOLATOR_HH_
#define PERFORMANCE_BENCHMARKS_BENCHMARK_INTERPOLATOR_HH_

#include "rascal/json_io.hh"
#include "rascal/math/interpolator.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "benchmarks.hh"

#include <functional>
#include <iostream>
#include <map>

namespace rascal {

  using internal::AtomicSmearingType;
  using internal::RadialBasisType;
  using internal::RadialContribution;

  using math::Hyp1f1;
  using math::InterpolatorMatrixUniformCubicSpline;
  using math::InterpolatorScalarUniformCubicSpline;
  using math::Matrix_Ref;
  using math::Matrix_t;
  using math::RefinementMethod_t;
  using math::Vector_Ref;
  using math::Vector_t;

  // For the random functionalities in the benchmarks
  static unsigned int SEED = 1597463007;  // 0x5f3759df

  /**
   * Computes the absolute error on the grid per grid point for scalar
   * interpolator
   */
  template <RefinementMethod_t RefinementMethod, class ErrorMethod>
  Vector_t compute_pointwise_absolute_grid_error(
      const std::shared_ptr<
          InterpolatorScalarUniformCubicSpline<RefinementMethod, ErrorMethod>> &
          intp) {
    Vector_t test_grid{intp->get_test_grid()};
    Vector_t test_grid_interpolated{intp->interpolate(Vector_Ref(test_grid))};
    Vector_t test_grid_evaluated{intp->eval(Vector_Ref(test_grid))};
    Vector_t error_grid{
        math::ErrorMethod<math::ErrorMetric_t::Absolute>::
            compute_pointwise_error(Vector_Ref(test_grid_interpolated),
                                    Vector_Ref(test_grid_evaluated))};
    return error_grid;
  }

  /**
   * Computes the absolute error on the grid per grid point for matrix
   * interpolator
   */
  template <RefinementMethod_t RefinementMethod, class ErrorMethod>
  Matrix_t compute_pointwise_absolute_grid_error(
      const std::shared_ptr<
          InterpolatorMatrixUniformCubicSpline<RefinementMethod, ErrorMethod>> &
          intp) {
    Vector_t test_grid{intp->get_test_grid()};
    Matrix_t test_grid_interpolated{
        intp->interpolate_to_vector(Vector_Ref(test_grid))};
    Matrix_t test_grid_evaluated{intp->eval(Vector_Ref(test_grid))};
    Matrix_t error_grid{
        math::ErrorMethod<math::ErrorMetric_t::Absolute>::
            compute_pointwise_error(Matrix_Ref(test_grid_interpolated),
                                    Matrix_Ref(test_grid_evaluated))};
    return error_grid;
  }

  class BaseInterpolatorDataset {
   public:
    enum class SupportedFunc { Identity, Gaussian, Hyp1f1 };
    enum class SupportedVecFunc { RadialContribution, Dummy };
  };

  /* -------------------- benchmark-parameters-start -------------------- */
  /**
   * Interpolator Dataset explanation
   *  static const json data = {
   *      // number of usages of the interpolation function on a grid, note this
   * is independent from the google benchmark parameter
   *      {"nbs_iterations", {1e3}},
   *
   *      // the range of the interpolator
   *      {"ranges", {std::make_pair(0, 16)}},
   *
   *      // the log function of the interpolator parameter error bound
   *      {"log_error_bounds", {-8}},
   *
   *      // the function used to interpolate
   *      {"func_names", {SupportedVecFunc::RadialContribution}},
   *
   *      // the parameters for `RadialContribution` or `SphericalExpansion` a
   *      // pair of (max_radial, max_angular)
   *      {"radial_angular", {std::make_pair(3, 3)},
   *
   *      // if the grid used for interpolation should be generated randomly,
   *      // this is performance critical when using the search method hunt or
   *      // locate in the interpolator
   *      {"random", {true}},
   *
   *      // parameter only for `SphericalExpansion`
   *      {"filenames",
   *       {"reference_data/CaCrP2O7_mvc-11955_symmetrized.json"}},
   *
   *      // parameter only for `SphericalExpansion` AdaptorStrict cutoff
   *      {"cutoffs", {8}}};
   */
  /* -------------------- benchmark-parameters-end -------------------- */

  /* -------------------- benchmark-dataset-start -------------------- */
  /**
   * Dataset for testing the spherical representations with interpolator. To
   * avoid code repetition `SphericalRepresentationBFixture` inherits from
   * `InterpolatorBFixture`. Because of the inheritance, all data parameters of
   * required by `InterpolatorBFixture` are also required by
   * `SphericalRepresentationBFixture`. But since not all of these parameters
   * are needed for the `SphericalRepresentationBFixture` we give them dummy
   * values.
   */
  struct SphericalDataset : public BaseInterpolatorDataset {
    using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
    static const json data() {
      static const json data = {
          {"nbs_iterations", {0}},             // dummy
          {"ranges", {std::make_pair(0, 0)}},  // dummy
          {"log_error_bounds", {-8, -10, -12}},
          {"func_names", {SupportedVecFunc::Dummy}},  // dummy
          {"radial_angular",
           {std::make_pair(3, 4), std::make_pair(6, 6), std::make_pair(8, 6)}},
          {"random", {true}},  // dummy
          {"filenames", {"../tests/reference_data/small_molecule.json"}},
          // please use only one file because google benchmark cant put strings
          // into their `Counter`, therefore the filename cannot be printed
          // {"filenames", {"../tests/reference_data/methane.json"}},
          {"cutoffs", {4}}};
      return data;
    }
  };
  /* -------------------- benchmark-dataset-end -------------------- */

  /**
   * RadialContributionDataset
   */
  struct RadialContributionDataset : public BaseInterpolatorDataset {
    using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
    static const json data() {
      static const json data = {
          {"nbs_iterations", {1e3}},
          {"ranges", {std::make_pair(0, 16)}},
          {"log_error_bounds", {-6, -8, -10, -12}},
          {"func_names", {SupportedVecFunc::RadialContribution}},
          {"radial_angular", {std::make_pair(6, 6)}},
          {"random", {true}}};
      return data;
    }
  };

  struct Hyp1f1Dataset : public BaseInterpolatorDataset {
    using SupportedFunc = typename BaseInterpolatorDataset::SupportedFunc;
    static const json data() {
      static const json data = {{"nbs_iterations", {1e3, 1e4, 1e5, 1e6}},
                                {"ranges", {std::make_pair(0, 16)}},
                                {"log_error_bounds", {-6, -8, -10, -12}},
                                {"func_names", {SupportedFunc::Hyp1f1}},
                                {"random", {true}}};
      return data;
    }
  };

  /**
   * Abstract class for all BFixtures using the interpolator to adapt for
   * different interpolator implementations.
   */
  template <class Dataset>
  class InterpolatorBFixture : public BaseBFixture<Dataset> {
   public:
    using Parent = BaseBFixture<Dataset>;
    using SupportedFunc = typename Dataset::SupportedFunc;

    void setup(const ::benchmark::State & state) {
      const json data = Dataset::data();
      // Because the two initialization processes share parameters of the
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
    void init_ref_points(const ::benchmark::State & state, const json & data) {
      auto range = this->template lookup<std::pair<double, double>>(
          data, "ranges", state);
      this->x1 = std::get<0>(range);
      this->x2 = std::get<1>(range);
      this->random = this->template lookup<bool>(data, "random", state);
      if (this->random) {
        math::Vector_t points_tmp =
            math::Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
        this->ref_points = math::Vector_t::Zero(this->nb_ref_points);
        for (int i{0}; i < this->ref_points.size(); i++) {
          this->ref_points(i) = points_tmp(rand_r(&SEED) % this->nb_ref_points);
        }
        for (int i{0}; i < this->ref_points.size(); i++) {
          assert(this->ref_points(i) >= this->x1 &&
                 this->ref_points(i) <= this->x2);
        }
      } else {
        this->ref_points =
            math::Vector_t::LinSpaced(this->nb_ref_points, this->x1, this->x2);
      }
    }
    virtual bool
    have_interpolator_parameters_changed(const ::benchmark::State & state,
                                         const json & data) const = 0;
    virtual void init_interpolator(const ::benchmark::State & state,
                                   const json & data) = 0;
  };

  /**
   * To cover different implementations of the interpolator, this abstract class
   * is the base class for all interpolator BFixtures.
   */
  template <class Dataset>
  class InterpolatorScalarBFixture : public InterpolatorBFixture<Dataset> {
   public:
    using Parent = InterpolatorBFixture<Dataset>;
    using SupportedFunc = typename Dataset::SupportedFunc;
    using Interpolator_t = math::InterpolatorScalarUniformCubicSpline<
        math::RefinementMethod_t::Exponential>;

    InterpolatorScalarBFixture<Dataset>() : Parent() {}

    std::shared_ptr<Interpolator_t> intp{};
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
      this->intp = std::make_shared<Interpolator_t>(
          this->func, this->x1, this->x2, this->error_bound);
    }

    void init_hyp1f1_function() {
      double n = 10;
      double l = 9;
      double a = 0.5 * (n + l + 3);
      double b = l + 1.5;
      // Important hyp1f1 a!=b, because for a==b the Hyp1f1 can be simplified
      // and computations does not take as long as it would on average
      auto hyp1f1 = math::Hyp1f1(a, b + 1, 200, 1e-15);
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
        throw std::runtime_error("SupportedFunc type is not known.");
        break;
      }
    }
  };

  /**
   * This class uses the matrix interpolator with the `RadialContribution`
   * function.
   */
  template <class Dataset>
  class InterpolatorMatrixBFixture : public InterpolatorBFixture<Dataset> {
   public:
    using Parent = InterpolatorBFixture<Dataset>;
    using SupportedVecFunc = typename Dataset::SupportedVecFunc;
    using Interpolator_t = math::InterpolatorMatrixUniformCubicSpline<
        math::RefinementMethod_t::Exponential>;

    InterpolatorMatrixBFixture<Dataset>() : Parent() {}

    std::shared_ptr<Interpolator_t> intp{};
    std::function<math::Matrix_t(double)> func{};
    SupportedVecFunc func_name{};
    int max_radial{0};
    int max_angular{0};
    std::shared_ptr<RadialContribution<RadialBasisType::GTO>> radial_contr{};

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

      auto new_radial_angular{this->template lookup<std::pair<int, int>>(
          data, "radial_angular", state)};
      int new_max_radial{std::get<0>(new_radial_angular)};
      int new_max_angular{std::get<1>(new_radial_angular)};
      return (have_scalar_parameters_changed ||
              new_max_radial != this->max_radial ||
              new_max_angular != this->max_angular);
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
      this->intp = std::make_shared<Interpolator_t>(
          this->func, this->x1, this->x2, this->error_bound, this->max_radial,
          this->max_angular + 1);
    }

    void init_radial_contribution_function(const ::benchmark::State & state,
                                           const json & data) {
      auto radial_angular = this->template lookup<std::pair<int, int>>(
          data, "radial_angular", state);
      this->max_radial = std::get<0>(radial_angular);
      this->max_angular = std::get<1>(radial_angular);
      json fc_hypers{{"type", "Constant"},
                     {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}};
      json hypers{
          {"gaussian_density", fc_hypers},
          {"max_radial", this->max_radial},
          {"max_angular", this->max_angular},
          {"cutoff_function", {{"cutoff", {{"value", 2.0}, {"unit", "A"}}}}}};
      // we cannot copy radial contribution to the lambda function, because
      // the copying has been disabled
      this->radial_contr =
          std::make_shared<RadialContribution<RadialBasisType::GTO>>(hypers);
      this->func = [&](double x) mutable {
        return this->radial_contr
            ->template compute_contribution<AtomicSmearingType::Constant>(x,
                                                                          0.5);
      };
    }

    void init_function(const ::benchmark::State & state, const json & data) {
      switch (this->func_name) {
      // This case uses the RadialContribution class as comparison and can
      // therefore directly be computed for different distances
      case SupportedVecFunc::RadialContribution:
        this->init_radial_contribution_function(state, data);
        break;
      default:
        throw std::runtime_error("SupportedVecFunc type is not known.");
        break;
      }
    }
  };

  /**
   * This class uses the RepresentationManager class as basis to call the
   * RadialContribution. Therefore benchmarks for different atomic structures
   * can be done.
   */
  template <class Representation_t, class Dataset>
  class SphericalRepresentationBFixture : public InterpolatorBFixture<Dataset> {
   public:
    // using Representation_t = CalculatorSphericalExpansion;
    using Manager_t = AdaptorStrict<AdaptorCenterContribution<
        AdaptorNeighbourList<StructureManagerCenters>>>;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Parent = InterpolatorBFixture<Dataset>;

    /**
     * @param use_interpolator determines if the SphericalExpansion should be
     * used with interpolator to estimate a time difference
     * @param compute_gradients determines if the SphericalExpansion should be
     * calculate gradients
     */
    SphericalRepresentationBFixture<Representation_t, Dataset>(
        bool use_interpolator, bool compute_gradients)
        : Parent(), use_interpolator{use_interpolator},
          compute_gradients{compute_gradients} {}

    bool use_interpolator;
    bool compute_gradients;
    int max_radial{0};
    int max_angular{0};
    std::string filename{};
    double cutoff{0};
    int nb_neighbours{0};
    ManagerPtr_t manager{};
    // to postpone initialization of representation from the time the object
    // is created, we create a list
    std::shared_ptr<Representation_t> representation_ptr{};
    json hypers{};

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
      // hyperparameters for the neighbourlist adaptor

      json ad1{{"name", "AdaptorNeighbourList"},
               {"initialization_arguments",
                {{"cutoff", this->cutoff},
                 {"consider_ghost_neighbours", false},
                 {"skin", 0.}}}};
      json ad2{{"name", "AdaptorCenterContribution"},
               {"initialization_arguments", {}}};
      // hyperparameters for the strict adaptor
      json ad3{{"name", "AdaptorStrict"},
               {"initialization_arguments", {{"cutoff", this->cutoff}}}};
      adaptors.emplace_back(ad1);
      adaptors.emplace_back(ad2);
      adaptors.emplace_back(ad3);

      AtomicStructure<3> atomic_structure{};
      atomic_structure.set_structure(this->filename);
      this->manager = make_structure_manager_stack<
          StructureManagerCenters, AdaptorNeighbourList,
          AdaptorCenterContribution, AdaptorStrict>(structure, adaptors);
      this->manager->update(atomic_structure);

      this->nb_neighbours = 0;
      for (auto center : this->manager) {
        for (auto _ : center) {
          this->nb_neighbours++;
        }
      }

      // make representation manager
      json hypers{{"max_radial", this->max_radial},
                  {"max_angular", this->max_angular},
                  {"soap_type", "PowerSpectrum"},
                  {"normalize", true},
                  {"compute_gradients", this->compute_gradients},
                  {"soap_type", "PowerSpectrum"}};

      json fc_hypers{{"type", "ShiftedCosine"},
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
      this->representation_ptr = std::make_shared<Representation_t>(hypers);
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
      double new_cutoff = this->template lookup<double>(data, "cutoffs", state);
      return (new_log_error_bound != this->log_error_bound ||
              new_max_radial != this->max_radial ||
              new_max_angular != this->max_angular ||
              new_filename != this->filename || new_cutoff != this->cutoff);
    }
  };

  template <class Dataset>
  using SphericalExpansionBFixture =
      SphericalRepresentationBFixture<CalculatorSphericalExpansion, Dataset>;
  template <class Dataset>
  using SphericalInvariantsBFixture =
      SphericalRepresentationBFixture<CalculatorSphericalInvariants, Dataset>;

}  // namespace rascal
#endif  // PERFORMANCE_BENCHMARKS_BENCHMARK_INTERPOLATOR_HH_
