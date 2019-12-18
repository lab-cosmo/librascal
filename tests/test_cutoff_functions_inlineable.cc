/**
 * @file   test_cutoff_functions_inlineable.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  Tests for the cutoff functions used in Behler-Parinello descriptors
 *
 * Copyright © 2019 Till Junge, Markus Stricker, LAMMM (EPFL), COSMO (EPFL)
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
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

#include "behler_fixtures.hh"
#include "test_structure.hh"

#include "rascal/representations/cutoff_functions_inlineable.hh"
#include "rascal/utils/units.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <exception>
#include <limits>

namespace rascal {
  using CutoffFunctions_t =
      boost::mpl::list<InlCutoffFunFixture<InlCutoffFunctionType::Cosine>>;

  BOOST_AUTO_TEST_SUITE(cutoff_functions_inlineable);

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(construction_test, Fix, CutoffFunctions_t,
                                   Fix) {
    using CutoffFun = typename Fix::CutoffFun;
    BOOST_CHECK_THROW(CutoffFun(this->unit_style, this->incorrect_put),
                      std::runtime_error);

    BOOST_CHECK_EQUAL(this->cut_fun.get_identifier(), "Cosine_1.1_Å");
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test, Fix, CutoffFunctions_t, Fix) {
    const double r_1{this->r_cut / 2};
    const double r_c{this->r_cut};

    BOOST_CHECK_EQUAL(this->cut_fun.f_c(r_c), 0);
    auto && grad{this->cut_fun.df_c(r_c)};

    // assuming that the ocutoff function is 1 at the centre, machine tol is the
    // correct comparison
    const auto tol{std::numeric_limits<double>::epsilon()};
    BOOST_CHECK_LT(std::abs(grad[0]), tol);
    BOOST_CHECK_LT(std::abs(grad[1]), tol);

    grad = this->cut_fun.df_c(r_1);
    BOOST_CHECK_GT(grad[0], 0);
    BOOST_CHECK_LT(grad[1], 0);

    BOOST_CHECK_EQUAL(grad[0], this->cut_fun.f_c(r_1));
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(compute_test, Fix, CutoffFunctions_t, Fix) {
    ManagerFixture<StructureManagerLammps> fixture{};
    auto strict{
        make_adapted_manager<AdaptorStrict>(fixture.manager, this->r_cut)};
    strict->update();
    constexpr bool ComputeDerivative{false};
    this->cut_fun.compute(*strict, ComputeDerivative);
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(compute_test_wrong_level, Fix,
                                   CutoffFunctions_t, Fix) {
    ManagerFixture<StructureManagerLammps> fixture{};
    auto strict{
        make_adapted_manager<AdaptorStrict>(fixture.manager, this->r_cut)};
    auto strict2{
        make_adapted_manager<AdaptorStrict>(strict, this->r_cut * .99)};
    strict2->update();
    constexpr bool ComputeDerivative{false};
    this->cut_fun.compute(*strict, ComputeDerivative);

    using Exception_t = std::runtime_error;
    BOOST_CHECK_THROW(this->cut_fun.compute(*strict2, ComputeDerivative),
                      Exception_t);
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(compute_derivative_test, Fix,
                                   CutoffFunctions_t, Fix) {
    ManagerFixture<StructureManagerLammps> fixture{};
    auto strict{
        make_adapted_manager<AdaptorStrict>(fixture.manager, this->r_cut)};
    strict->update();
    constexpr bool ComputeDerivative{true};
    this->cut_fun.compute(*strict, ComputeDerivative);
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
