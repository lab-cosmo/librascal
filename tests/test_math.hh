/**
 * @file   test_math.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief Test implementation of math functions
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef TESTS_TEST_MATH_HH_
#define TESTS_TEST_MATH_HH_

#include "json_io.hh"
#include "math/bessel.hh"
#include "math/gauss_legendre.hh"
#include "math/hyp1f1.hh"
#include "math/math_utils.hh"
#include "math/spherical_harmonics.hh"
#include "rascal_utility.hh"

#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <string>

namespace rascal {

  struct SphericalHarmonicsRefFixture {
    SphericalHarmonicsRefFixture() {
      json ref_data;
      std::ifstream ref_file(this->ref_filename);
      ref_file >> ref_data;
      unit_vectors = ref_data.at("unit_vectors").get<StdVector2Dim_t>();
      harmonics = ref_data.at("harmonics").get<StdVector3Dim_t>();
      // a vector which contains Switzerland's mountain landscape
      alps = ref_data.at("alps").get<StdVector3Dim_t>();
    }

    ~SphericalHarmonicsRefFixture() = default;

    std::string ref_filename = "reference_data/spherical_harmonics_test.json";

    using StdVector2Dim_t = std::vector<std::vector<double>>;
    using StdVector3Dim_t = std::vector<std::vector<std::vector<double>>>;
    StdVector2Dim_t unit_vectors{};
    StdVector3Dim_t harmonics{};
    StdVector3Dim_t alps{};
    bool verbose{false};
  };

  struct SphericalHarmonicsClassRefFixture {
    SphericalHarmonicsClassRefFixture() {
      this->ref_data =
          json::from_ubjson(internal::read_binary_file(this->ref_filename));
    }

    ~SphericalHarmonicsClassRefFixture() = default;

    std::string ref_filename =
        "reference_data/spherical_harmonics_reference.ubjson";
    json ref_data{};
    // TODO(alex) replace this with one variable VerbosityValues verbosity
    // for general test information
    bool info{false};
    // for detailed tests information of computed values
    bool verbose{false};
  };

  struct GaussLegendreRefFixture {
    GaussLegendreRefFixture() {
      this->ref_data =
          json::from_ubjson(internal::read_binary_file(this->ref_filename));
    }

    ~GaussLegendreRefFixture() = default;

    std::string ref_filename = "reference_data/gauss_legendre_reference.ubjson";

    json ref_data{};
    bool verbose{false};
  };

  struct ModifiedBesselFirstKindRefFixture {
    ModifiedBesselFirstKindRefFixture() {
      this->ref_data = json_io::load_txt(this->ref_filename);
    }

    ~ModifiedBesselFirstKindRefFixture() = default;

    std::string ref_filename =
        "reference_data/modified_bessel_first_kind_reference.json";

    json ref_data{};
    math::ModifiedSphericalBessel j_v_complete_square{};
    bool verbose{true};
  };

  /**
   * Wrapper of the SphericalHarmonics class to interface with the gradient
   * tester
   *
   * See the documentation for test_gradients() below; the calculator object
   * passed to it must provide the functions f() and grad_f() as below.
   */
  template <size_t max_angular>
  struct SphericalHarmonicsGradientsCalculator {
    SphericalHarmonicsGradientsCalculator()
        : harmonics_calculator{math::SphericalHarmonics(true)} {}

    static const size_t n_arguments = 3;
    using Matrix_Ref = typename Eigen::Ref<const math::Matrix_t>;
    using Vector_Ref = typename Eigen::Ref<const math::Vector_t>;

    void precompute() { this->harmonics_calculator.precompute(max_angular); }

    Vector_Ref f(const Eigen::Vector3d & inputs_v) {
      Eigen::Vector3d my_inputs = inputs_v / inputs_v.norm();
      this->harmonics_calculator.calc(my_inputs, false);
      return this->harmonics_calculator.get_harmonics();
    }
    Matrix_Ref grad_f(const Eigen::Vector3d & inputs_v) {
      this->harmonics_calculator.calc(inputs_v, true);
      return this->harmonics_calculator.get_harmonics_derivatives();
    }

    math::SphericalHarmonics harmonics_calculator;
  };

  /**
   * Fixture for testing a the gradient of a real function of N real arguments
   *
   * (Verifies that the gradient is the same as the converged value of the
   * finite-difference approximation along the given directions)
   *
   * Parameters should be provided in a JSON input file, as follows:
   *
   * @param function_inputs List of vectors of function arguments at which to
   *                        test the gradient
   *
   * @param direction_mode How the finite-difference directions are specified;
   *                       options are "Cartesian" (once along each independent
   *                       argument), "Random" (exactly what it says on the
   *                       tin), and "Provided" (given in input file, see below)
   *
   * @param displacement_directions List of vectors along which to displace the
   *                                inputs, in case "direction_mode" is
   *                                "Provided".
   */
  class GradientTestFixture {
   public:
    explicit GradientTestFixture(std::string input_filename) {
      using Eigen::ArrayXd;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      json input_data;

      std::ifstream input_file(input_filename);
      input_file >> input_data;
      this->function_inputs =
          input_data.at("function_inputs").get<StdVector2Dim_t>();
      this->n_arguments = function_inputs[0].size();

      this->displacement_directions =
          this->get_displacement_directions(input_data, this->n_arguments);
      this->verbosity = get_verbosity(input_data);
      if (input_data.find("fd_error_tol") != input_data.end()) {
        this->fd_error_tol = input_data["fd_error_tol"].get<double>();
      }
    }

    ~GradientTestFixture() = default;

    static Eigen::MatrixXd get_displacement_directions(json & input_data,
                                                       size_t n_arguments) {
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      MatrixXd directions;
      std::string direction_mode =
          input_data.at("direction_mode").get<std::string>();
      if (direction_mode.compare("Cartesian") == 0) {
        directions = MatrixXd::Identity(n_arguments, n_arguments);
      } else if (direction_mode.compare("Random") == 0) {
        size_t n_directions{input_data.at("n_directions").get<size_t>()};
        directions = MatrixXd::Random(n_directions, n_arguments);
        directions.rowwise().normalize();
      } else if (direction_mode.compare("Provided") == 0) {
        StdVector2Dim_t directions_in =
            input_data.at("displacement_directions").get<StdVector2Dim_t>();
        directions.resize(directions_in.size(), n_arguments);
        int row_idx{0};
        for (auto & direction : directions_in) {
          directions.row(row_idx++) =
              Eigen::Map<VectorXd>(direction.data(), n_arguments);
        }
        directions.rowwise().normalize();
      } else {
        std::cerr << "Unknown direction mode \'" << direction_mode;
        std::cerr << "\', assuming Cartesian" << std::endl;
        directions = MatrixXd::Identity(n_arguments, n_arguments);
      }
      return directions;
    }

    enum struct VerbosityValue {
      NORMAL = 0,  // Print nothing
      INFO = 10,   // Print one line of info for each gradient step
      DEBUG = 20   // Print as much as possible
    };

    static VerbosityValue get_verbosity(json & input_data) {
      VerbosityValue verbosity_in{VerbosityValue::NORMAL};
      if (input_data.find("verbosity") != input_data.end()) {
        std::string verbosity_str = input_data["verbosity"].get<std::string>();
        if (verbosity_str.compare("INFO") == 0) {
          verbosity_in = VerbosityValue::INFO;
        } else if (verbosity_str.compare("DEBUG") == 0) {
          verbosity_in = VerbosityValue::DEBUG;
        } else if (verbosity_str.compare("NORMAL") == 0) {
          verbosity_in = VerbosityValue::NORMAL;
        } else {
          std::cerr << "Unknown verbosity value \'" << verbosity_str;
          std::cerr << "\', assuming NORMAL" << std::endl;
          verbosity_in = VerbosityValue::NORMAL;
        }
      }
      return verbosity_in;
    }

    using StdVector2Dim_t = std::vector<std::vector<double>>;

    StdVector2Dim_t function_inputs{};
    Eigen::MatrixXd displacement_directions{};
    size_t n_arguments{0};

    /**
     * The error of the finite-difference against the analytical derivatives
     * isn't going to be arbitrarily small, due to the interaction of
     * finite-difference and finite precision effects.  The default here
     * reflects this, while allowing different tests to change it -- keeping in
     * mind that the automated test is really more intended to be a sanity
     * check on the implementation than a rigorous convergence test.
     */
    double fd_error_tol{1E-6};
    VerbosityValue verbosity{VerbosityValue::NORMAL};

   protected:
    GradientTestFixture() {}
  };

  namespace internal {
    /**
     * Template to select the Eigen type for the argument vector: Fixed or
     * dynamic size?  Fixed-size template below, dynamic-size default here.
     */
    template <typename Calculator_t, typename = void>
    struct Argument_t {
      typedef Eigen::VectorXd type_vec;
      typedef Eigen::MatrixXd type_mat;
    };

    // Works by SFINAE -- if Calculator has no static integral member
    // `n_arguments`, this specialization fails and the default (above) is
    // selected
    template <typename Calculator_t>
    struct Argument_t<Calculator_t, std::enable_if_t<std::is_integral<decltype(
                                        Calculator_t::n_arguments)>::value>> {
      typedef Eigen::Matrix<double, Calculator_t::n_arguments, 1> type_vec;
      typedef Eigen::Matrix<double, Calculator_t::n_arguments, Eigen::Dynamic>
          type_mat;
    };
  }  // namespace internal

  /**
   * Numerically verify that a given function and its gradient are consistent
   *
   * @param function_calculator An object that provides both the function and
   *                            its gradient
   *
   * @param params Test fixture object, e.g a GradientTestFixture or something
   *               providing the same information (i.e. function_inputs,
   *               displacement_directions, n_arguments, and verbosity)
   *
   * The function_calculator object may be of any type, as long as it provides
   * two functions, f() and grad_f(), to calculate the function and its gradient
   * (derivative for functions with one input, Jacobian for functions with
   * multiple outputs -- the output dimension is expected to correspond to
   * columns).  Both functions must accept an Eigen::Vector, corresponding to
   * the function input, of dimension determined in the data file (read by
   * GradientTestFixture).  This function additionally guarantees that f() will
   * be called before grad_f() with the same input.
   *
   * If the functions f() and grad_f() are designed to accept fixed-size vectors
   * (i.e. if the size of the argument is known at compile time), be sure to
   * define, in the FunctionProvider class, a
   * `constexpr static size_t n_arguments` member with the size of the argument
   * vector.  This will ensure that the gradient tester uses the corresponding
   * fixed-size Eigen vectors/matrices as inputs.
   */
  template <typename FunctionProvider_t, typename TestFixture_t>
  void test_gradients(FunctionProvider_t function_calculator,
                      TestFixture_t params) {
    using ArgumentVec_t =
        typename internal::Argument_t<FunctionProvider_t>::type_vec;
    using ArgumentMat_t =
        typename internal::Argument_t<FunctionProvider_t>::type_mat;

    using VerbosityValue = typename GradientTestFixture::VerbosityValue;

    for (auto inputs : params.function_inputs) {
      ArgumentVec_t argument_vector =
          Eigen::Map<ArgumentVec_t>(inputs.data(), params.n_arguments, 1);
      Eigen::MatrixXd values = function_calculator.f(argument_vector);
      ArgumentMat_t jacobian = function_calculator.grad_f(argument_vector);
      if (params.verbosity >= VerbosityValue::INFO) {
        std::cout << std::string(30, '-') << std::endl;
        std::cout << "Input vector: " << argument_vector.transpose();
        std::cout << std::endl;
      }
      if (params.verbosity >= VerbosityValue::DEBUG) {
        std::cout << "Function values:" << values << std::endl;
        std::cout << "Jacobian:" << jacobian << std::endl;
      }
      for (int disp_idx{0}; disp_idx < params.displacement_directions.rows();
           disp_idx++) {
        ArgumentVec_t displacement_direction =
            params.displacement_directions.row(disp_idx);
        // Compute the directional derivative(s)
        Eigen::VectorXd directional =
            displacement_direction.transpose() * jacobian;
        if (params.verbosity >= VerbosityValue::INFO) {
          std::cout << "FD direction: " << displacement_direction.transpose();
          std::cout << std::endl;
        }
        if (params.verbosity >= VerbosityValue::DEBUG) {
          std::cout << "Analytical derivative: " << directional.transpose();
          std::cout << std::endl;
        }
        // TODO(max) this inner loop is a good candidate to move to its own
        // function... if we can pass along those argument typedefs as well
        double min_error{HUGE_VAL};
        for (double dx = 1E-2; dx > 1E-10; dx *= 0.1) {
          if (params.verbosity >= VerbosityValue::INFO) {
            std::cout << "dx = " << dx << "\t";
          }
          ArgumentVec_t displacement = dx * displacement_direction;
          // Compute the finite-difference derivative using a
          // centred-difference approach
          Eigen::VectorXd fun_plus{
              function_calculator.f(argument_vector + displacement)};
          Eigen::VectorXd fun_minus{
              function_calculator.f(argument_vector - displacement)};
          Eigen::VectorXd fd_derivatives = 0.5 / dx * (fun_plus - fun_minus);
          double fd_error{0.};
          double fd_quotient{0.};
          size_t nonzero_count{0};
          for (int dim_idx{0}; dim_idx < fd_derivatives.size(); dim_idx++) {
            if (std::abs(directional(dim_idx)) < 10 * math::dbl_ftol) {
              fd_error += fd_derivatives(dim_idx);
            } else {
              fd_quotient += (fd_derivatives(dim_idx) / directional(dim_idx));
              fd_error += (fd_derivatives(dim_idx) - directional(dim_idx)) /
                          directional(dim_idx);
              ++nonzero_count;
            }
          }
          if (nonzero_count > 0) {
            fd_quotient = fd_quotient / nonzero_count;
          }
          fd_error = fd_error / fd_derivatives.size();
          if (params.verbosity >= VerbosityValue::INFO) {
            std::cout << "Average rel FD error: " << fd_error << "\t";
            std::cout << "Average FD quotient:  " << fd_quotient << std::endl;
          }
          if (std::abs(fd_error) < min_error) {
            min_error = std::abs(fd_error);
          }
          if (params.verbosity >= VerbosityValue::DEBUG) {
            Eigen::VectorXd fd_error_cwise = (fd_derivatives - directional);
            std::cout << "error            = " << fd_error_cwise.transpose();
            std::cout << std::endl;
            std::cout << "(FD derivative   = " << fd_derivatives.transpose();
            std::cout << ")" << std::endl;
          }
        }  // for (double dx...) (displacement magnitudes)
        BOOST_CHECK_SMALL(min_error, params.fd_error_tol);
      }  // for (int disp_idx...) (displacement directions)
    }    // for (auto inputs...) (function inputs)
  }

  struct Hyp1F1RefFixture {
    Hyp1F1RefFixture() {
      this->ref_data =
          json::from_ubjson(internal::read_binary_file(this->ref_filename));
    }

    ~Hyp1F1RefFixture() = default;

    std::string ref_filename = "reference_data/hyp1f1_reference.ubjson";

    json ref_data{};
    bool verbose{false};
  };

  struct Hyp1f1SphericalExpansionFixture {
    Hyp1f1SphericalExpansionFixture() {
      for (auto & l_max : l_maxs) {
        for (auto & n_max : n_maxs) {
          hyp1f1.emplace_back(false, 1e-14);
          hyp1f1.back().precompute(n_max, l_max);
          hyp1f1_recursion.emplace_back(true, 1e-14);
          hyp1f1_recursion.back().precompute(n_max, l_max);
        }
      }

      for (auto & rc : rcs) {
        facs_b.emplace_back();
        for (size_t il{0}; il < l_maxs.size();) {
          for (auto & n_max : n_maxs) {
            facs_b.back().emplace_back(n_max);
            for (int n{0}; n < n_max; ++n) {
              double sigma_n{(rc - smooth_width) * std::max(std::sqrt(n), 1.) /
                             n_max};
              facs_b.back().back()(n) = 0.5 * math::pow(sigma_n, 2);
            }
          }
          il++;
        }
      }
    }

    ~Hyp1f1SphericalExpansionFixture() = default;

    std::vector<int> l_maxs{{4, 5, 9, 15, 16, 20}};
    std::vector<int> n_maxs{{4, 5, 9, 15, 16, 20}};
    std::vector<math::Hyp1f1SphericalExpansion> hyp1f1{};
    std::vector<math::Hyp1f1SphericalExpansion> hyp1f1_recursion{};
    std::vector<std::vector<Eigen::VectorXd>> facs_b{};
    std::vector<double> r_ijs{1., 2., 3., 4., 5.5, 6.5, 7.5, 7.9};
    std::vector<double> fac_as{0.4};
    std::vector<double> rcs{2., 3., 5., 7., 8.};
    double smooth_width{0.5};
    bool verbose{false};
  };

  struct Hyp1f1GradientProvider {
    Hyp1f1GradientProvider(size_t max_radial, size_t max_angular, double fac_a,
                           Eigen::Ref<Eigen::VectorXd> fac_b)
        : max_radial{max_radial}, max_angular{max_angular}, fac_a{fac_a} {
      this->fac_b.resize(max_angular, 1);
      this->fac_b = fac_b;
      this->hyp1f1_calculator.precompute(max_radial, max_angular);
    }

    ~Hyp1f1GradientProvider() = default;

    Eigen::Ref<Eigen::Array<double, 1, Eigen::Dynamic>>
    f(const Eigen::Matrix<double, 1, 1> & input_v) {
      this->hyp1f1_calculator.calc(input_v(0), this->fac_a, this->fac_b);
      Eigen::MatrixXd result(this->max_radial, this->max_angular + 1);
      result = this->hyp1f1_calculator.get_values();
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    Eigen::Ref<Eigen::Array<double, 1, Eigen::Dynamic>>
    grad_f(const Eigen::Matrix<double, 1, 1> & input_v) {
      this->hyp1f1_calculator.calc(input_v(0), this->fac_a, this->fac_b, true);
      Eigen::MatrixXd result(this->max_radial, this->max_angular + 1);
      result = this->hyp1f1_calculator.get_derivatives();
      Eigen::Map<Eigen::Array<double, 1, Eigen::Dynamic>> result_flat(
          result.data(), 1, result.size());
      return result_flat;
    }

    math::Hyp1f1SphericalExpansion hyp1f1_calculator{true, 1e-13, 200};
    size_t max_radial;
    size_t max_angular;
    double fac_a{};
    Eigen::VectorXd fac_b{};
  };

  template <class CutoffFunction>
  struct CutoffGradientProvider {
    explicit CutoffGradientProvider(CutoffFunction & cutoff)
        : cutoff_calculator{cutoff} {}

    ~CutoffGradientProvider() = default;

    Eigen::MatrixXd f(const Eigen::Matrix<double, 1, 1> & input_v) {
      Eigen::MatrixXd result(1, 1);
      result(0) = this->cutoff_calculator.f_c(input_v(0));
      return result;
    }

    Eigen::MatrixXd grad_f(const Eigen::Matrix<double, 1, 1> & input_v) {
      Eigen::MatrixXd result(1, 1);
      result(0) = this->cutoff_calculator.df_c(input_v(0));
      return result;
    }
    static const size_t n_arguments = 1;
    CutoffFunction & cutoff_calculator;
  };

}  // namespace rascal

#endif  // TESTS_TEST_MATH_HH_
