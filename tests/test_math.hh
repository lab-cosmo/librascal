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
    GradientTestFixture() {
      using Eigen::ArrayXd;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      json input_data;

      std::ifstream input_file(this->input_filename);
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
    }

    ~GradientTestFixture() = default;

    using StdVector2Dim_t = std::vector<std::vector<double>>;
    StdVector2Dim_t function_inputs{};
    std::vector<double> displacement_lengths{};
    Eigen::MatrixXd displacement_directions{};
    std::string input_filename{
      "reference_data/spherical_harmonics_gradient_test.json"};
    size_t n_arguments{0};
    bool verbose{true};
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
}  // namespace rascal

#endif  // TESTS_TEST_MATH_HH_
