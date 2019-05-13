/**
 * file   test_math.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
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

#include "tests.hh"
#include "json_io.hh"
#include "math/math_interface.hh"
#include "math/math_utils.hh"

#include <fstream>
#include <string>
#include <Eigen/Dense>

namespace rascal {

  struct ManagerFixtureMath {
    ManagerFixtureMath()
        : numbers(4, 3), results_hyp2f1(3), results_airy(3, 4) {
      // clang-format off
      numbers <<   1, 0.1,   2,
                   1,   3,   9,
                   2,   7,   6,
                 0.5, 0.2, 0.3;

      results_hyp2f1 << 1.3862943611198901,
                        1.0090833356005495,
                        3.0875740550280937;
      results_airy <<
        0.13529241631288147, -0.15914744129679328,  1.2074235949528715,
         0.9324359333927756,    0.329203129943538, -0.2571304219075862,
          0.659861690194189,  0.45151263114964657,  0.03492413042327436,
       -0.05309038443365388,   3.2980949999782143,  4.10068204993289;
      // clang-format on
    }

    ~ManagerFixtureMath() {}

    Eigen::Matrix<double, 4, Eigen::Dynamic> numbers;
    Eigen::Matrix<double, 1, Eigen::Dynamic> results_hyp2f1;
    Eigen::Matrix<double, 3, Eigen::Dynamic> results_airy;
    bool vebose{false};
  };

  struct SphericalHarmonicsRefFixture {
    SphericalHarmonicsRefFixture() {
      json ref_data;
      std::ifstream ref_file(this->ref_filename);
      ref_file >> ref_data;
      unit_vectors = ref_data.at("unit_vectors").get<StdVector2Dim_t>();
      harmonics = ref_data.at("harmonics").get<StdVector3Dim_t>();
      alps = ref_data.at("alps").get<StdVector3Dim_t>();
    }

    ~SphericalHarmonicsRefFixture() = default;

    std::string ref_filename = "reference_data/spherical_harmonics_test.json";

    using StdVector2Dim_t = std::vector<std::vector<double>>;
    using StdVector3Dim_t = std::vector<std::vector<std::vector<double>>>;
    StdVector2Dim_t unit_vectors{};
    StdVector3Dim_t harmonics{};
    StdVector3Dim_t alps{};
    bool verbose{true};
  };

  /**
   * Fixture for testing a the gradient of a scalar function of N real arguments
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
  struct GradientTestFixture {
    explicit GradientTestFixture(std::string input_filename) {
      using Eigen::ArrayXd;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      json input_data;

      std::ifstream input_file(input_filename);
      input_file >> input_data;
      function_inputs = input_data.at("function_inputs").get<StdVector2Dim_t>();
      n_arguments = function_inputs[0].size();

      // Get/make the displacement directions
      std::string direction_mode = input_data.at(
          "direction_mode").get<std::string>();
      if (direction_mode.compare("Cartesian") == 0) {
        displacement_directions = MatrixXd::Identity(n_arguments, n_arguments);
      } else if (direction_mode.compare("Random") == 0) {
        size_t n_directions{input_data.at("n_directions").get<size_t>()};
        displacement_directions = MatrixXd::Random(n_directions, n_arguments);
        displacement_directions.rowwise().normalize();
      } else if (direction_mode.compare("Provided") == 0) {
        StdVector2Dim_t directions_in = input_data.at(
            "displacement_directions").get<StdVector2Dim_t>();
        displacement_directions.resize(directions_in.size(), n_arguments);
        int row_idx{0};
        for (auto it{directions_in.begin()}; it != directions_in.end(); it++) {
          displacement_directions.row(row_idx++) =
            Eigen::Map<VectorXd> (it->data(), 1, n_arguments);
        }
        displacement_directions.rowwise().normalize();
      } else {
        std::cerr << "Unknown direction mode \'" << direction_mode;
        std::cerr << "\', assuming Cartesian" << std::endl;
        displacement_directions = MatrixXd::Identity(n_arguments, n_arguments);
      }

      if (input_data.find("verbose") != input_data.end()) {
        verbose = input_data["verbose"].get<bool>();
      }
    }

    ~GradientTestFixture() = default;

    using StdVector2Dim_t = std::vector<std::vector<double>>;
    StdVector2Dim_t function_inputs{};
    std::vector<double> displacement_lengths{};
    Eigen::MatrixXd displacement_directions{};
    //std::string input_filename{
      //"reference_data/spherical_harmonics_gradient_test.json"};
    size_t n_arguments{0};
    bool verbose{false};
  };

  template<int max_angular>
  class SphericalHarmonicsWithGradients {
   public:
    SphericalHarmonicsWithGradients() = default;

    ~SphericalHarmonicsWithGradients() = default;

    Eigen::Array<double, 1, (max_angular+1)*(max_angular+1)>
    f(const Eigen::Vector3d & inputs_v) {
      // Renormalize the inputs to project out the r gradients
      Eigen::Vector3d my_inputs = inputs_v / inputs_v.norm();
      return math::compute_spherical_harmonics(my_inputs, max_angular);
    }

    Eigen::Array<double, 3, (max_angular+1)*(max_angular+1)>
    grad_f(const Eigen::Vector3d & inputs_v) {
      Eigen::Array<double, 4, (max_angular+1)*(max_angular+1)>
        harmonics_derivatives{math::compute_spherical_harmonics_derivatives(
                                                        inputs_v, max_angular)};
      // The gradients are with respect to the central atom, but the
      // displacements here affect the neighbouring atom, hence the sign
      return -1.0 * harmonics_derivatives.bottomRows(3);
    }
  };

  template<typename FunctionProvider_t> void test_gradients(
                                              std::string data_filename) {
    FunctionProvider_t function_calculator{};
    GradientTestFixture params{data_filename};
    Eigen::MatrixXd values;
    Eigen::MatrixXd jacobian;
    Eigen::RowVectorXd argument_vector;
    Eigen::VectorXd displacement_direction;
    Eigen::Vector3d displacement;
    Eigen::MatrixXd directional;
    Eigen::MatrixXd fd_derivatives;
    Eigen::MatrixXd fd_error_cwise;
    // This error isn't going to be arbitrarily small, due to the interaction of
    // finite-difference and finite precision effects.  Just set it to something
    // reasonable and check it explicitly if you really want to be sure (paying
    // attention to the change of the finite-difference gradients from one step
    // to the next).  The automated test is really more intended to be a sanity
    // check on the implementation anyway.
    constexpr double fd_error_tol = 1E-8;
    for (auto inputs_it{params.function_inputs.begin()};
         inputs_it != params.function_inputs.end(); inputs_it++) {
      argument_vector = Eigen::Map<Eigen::RowVectorXd>
                                    (inputs_it->data(), 1, params.n_arguments);
      values = function_calculator.f(argument_vector);
      jacobian = function_calculator.grad_f(argument_vector);
      std::cout << std::string(30, '-') << std::endl;
      std::cout << "Direction vector: " << argument_vector << std::endl;
      if (params.verbose) {
        std::cout << "Values:" << values << std::endl;
        std::cout << "Jacobian:" << jacobian << std::endl;
      }
      for (int disp_idx{0}; disp_idx < params.displacement_directions.rows();
           disp_idx++) {
        displacement_direction = params.displacement_directions.row(disp_idx);
        // Compute the directional derivative(s)
        directional = displacement_direction.adjoint() * jacobian;
        std::cout << "FD direction: " << displacement_direction.adjoint();
        std::cout << std::endl;
        if (params.verbose) {
          std::cout << "Analytical derivative: " << directional << std::endl;
        }
        double min_error{HUGE_VAL};
        Eigen::MatrixXd fd_last{Eigen::MatrixXd::Zero(1, directional.size())};
        for (double dx = 1E-2; dx > 1E-10; dx *= 0.1) {
          std::cout << "dx = " << dx << "\t";
          displacement = dx * displacement_direction;
          // Compute the finite-difference derivative using a
          // centred-difference approach
          fd_derivatives = 0.5 / dx * (
            function_calculator.f(argument_vector + displacement.adjoint())
          - function_calculator.f(argument_vector - displacement.adjoint()));
          double fd_error{0.};
          double fd_quotient{0.};
          size_t nonzero_count{0};
          for (int dim_idx{0}; dim_idx < fd_derivatives.size(); dim_idx++) {
            if (std::abs(directional(dim_idx)) < 10*math::dbl_ftol) {
              fd_error += fd_derivatives(dim_idx);
            } else {
              fd_quotient += (fd_derivatives(dim_idx) / directional(dim_idx));
              fd_error += (fd_derivatives(dim_idx) - directional(dim_idx))
                          / directional(dim_idx);
              ++nonzero_count;
            }
          }
          if (nonzero_count > 0) {
            fd_quotient = fd_quotient / nonzero_count;
          }
          fd_error = fd_error / fd_derivatives.size();
          std::cout << "Average rel FD error: " << fd_error << "\t";
          std::cout << "Average FD quotient:  " << fd_quotient << std::endl;
          if (std::abs(fd_error) < min_error) {min_error = std::abs(fd_error);}
          if (params.verbose) {
            fd_error_cwise = (fd_derivatives - directional);
            std::cout << "error            = " << fd_error_cwise << std::endl;
            std::cout << "(FD derivative   = " << fd_derivatives << ")";
            std::cout << std::endl;
            std::cout << "(minus last step = " << fd_derivatives - fd_last;
            std::cout << ")" << std::endl;
          }
          fd_last = fd_derivatives;
        } // for (double dx...) (displacement magnitudes)
        BOOST_CHECK_SMALL(min_error, fd_error_tol);
      } // for (int disp_idx...) (displacement directions)
    } // for (auto inputs_it...) (function inputs)
  }
}  // namespace rascal

#endif  // TESTS_TEST_MATH_HH_
