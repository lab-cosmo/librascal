/**
 * @file   test_cutoff_functions.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   21 May 2019
 *
 * @brief Test the implementation of cutoff functions
 *
 * Copyright © 2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "test_math.hh"

#include "rascal/representations/cutoff_functions.hh"

#include <boost/test/unit_test.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(cutoff_tests);

  /* ---------------------------------------------------------------------- */
  /**
   * Test the gradients of the shifted cosine cutoff function
   */
  BOOST_AUTO_TEST_CASE(shifted_cosine_gradient_test) {
    static const bool verbose{false};
    std::vector<json> fc_hypers{{{"cutoff_function_type", "ShiftedCosine"},
                                 {"interaction_cutoff", 3.0},
                                 {"cutoff_smooth_width", 0.5}}};
    using Cutoff_t =
        internal::CutoffFunction<internal::CutoffFunctionType::ShiftedCosine>;

    for (auto & fc_hyper : fc_hypers) {
      if (verbose) {
        std::cout << fc_hyper.dump(2) << std::endl;
      }

      Cutoff_t cutoff(fc_hyper);

      CutoffGradientProvider<Cutoff_t> cutoff_calculator(cutoff);

      GradientTestFixture fix{
          "reference_data/tests_only/cutoff_function_test.json"};
      test_gradients(cutoff_calculator, fix);
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test the gradients of the radial scaling cutoff function
   */
  BOOST_AUTO_TEST_CASE(radial_scaling_gradient_test) {
    static const bool verbose{false};
    std::vector<json> fc_hypers{
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 1}, {"scale", 2}, {"exponent", 3}}}},
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 0}, {"scale", 2}, {"exponent", 3}}}},
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 0}, {"scale", 2}, {"exponent", 0}}}}};

    using Cutoff_t =
        internal::CutoffFunction<internal::CutoffFunctionType::RadialScaling>;

    for (auto & fc_hyper : fc_hypers) {
      if (verbose) {
        std::cout << fc_hyper.dump(2) << std::endl;
      }

      Cutoff_t cutoff(fc_hyper);

      CutoffGradientProvider<Cutoff_t> cutoff_calculator(cutoff);

      GradientTestFixture fix{
          "reference_data/tests_only/cutoff_function_test.json"};
      test_gradients(cutoff_calculator, fix);
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test if the exponent is zero case (m = 0) is correctly identified in the
   * radial scaling cutoff function
   */
  BOOST_AUTO_TEST_CASE(radial_scaling_exponent_zero_test) {
    static const bool verbose{false};
    std::vector<json> fc_hypers{
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 1}, {"scale", 1}, {"exponent", 0}}}},
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 0}, {"scale", 0}, {"exponent", 0}}}},
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 0}, {"scale", 1}, {"exponent", 0}}}},
        {{"cutoff_function_type", "RadialScaling"},
         {"interaction_cutoff", 3.0},
         {"cutoff_smooth_width", 0.5},
         {"cutoff_function_parameters",
          {{"rate", 1}, {"scale", 0}, {"exponent", 0}}}}};
    using Cutoff_t =
        internal::CutoffFunction<internal::CutoffFunctionType::RadialScaling>;

    for (auto & fc_hyper : fc_hypers) {
      if (verbose) {
        std::cout << fc_hyper.dump(2) << std::endl;
      }

      Cutoff_t cutoff(fc_hyper);
      BOOST_CHECK_EQUAL(cutoff.value(2), 1);
      BOOST_CHECK_EQUAL(cutoff.grad(2), 0);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
