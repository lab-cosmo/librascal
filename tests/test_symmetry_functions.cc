/**
 * @file   test_symmetry_functions.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  Tests for the symmetry functions used in Behler-Parinello descriptors
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

#include "rascal/representations/symmetry_functions.hh"
#include "behler_fixtures.hh"
#include "rascal/utils/json_io.hh"
#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <exception>
#include <limits>

#include "Eigen/Dense"

namespace rascal {
  using SymmetryFunctions_t =
      boost::mpl::list<SymmetryFunFixture<SymmetryFunctionType::Gaussian>>;

  BOOST_AUTO_TEST_SUITE(symmetry_functions_behler);

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(construction_test, Fix, SymmetryFunctions_t,
                                   Fix) {
    using SymFun = typename Fix::SymFun;

    BOOST_CHECK_THROW(SymFun(this->unit_style, this->incorrect_put),
                      std::runtime_error);
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test, Fix, SymmetryFunctions_t, Fix) {
    using Vec_t = Eigen::Matrix<double, ThreeD, 1>;
    const Vec_t n_ij{[]() {
      Vec_t vec{Vec_t::Random()};
      return vec / vec.norm();
    }()};

    BOOST_CHECK_NO_THROW(this->sym_fun.f_sym(this->r_ij));
    BOOST_CHECK_NO_THROW(this->sym_fun.df_sym(this->r_ij, n_ij));

    double f1{this->sym_fun.f_sym(this->r_ij)};
    double f2{std::get<0>(this->sym_fun.df_sym(this->r_ij, n_ij))};

    BOOST_CHECK_EQUAL(f1, f2);
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
