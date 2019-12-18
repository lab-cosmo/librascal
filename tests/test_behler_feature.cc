/**
 * @file   test_behler_feature.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  testing Behler-Parinello G-functions
 *
 * Copyright © 2019 Till Junge
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

// #include "rascal/representations/behler_feature.hh"
#include "behler_fixtures.hh"
#include "test_structure.hh"

#include "rascal/utils/json_io.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {

  template <SymmetryFunctionType SymFunType, InlCutoffFunctionType CutFunType>
  struct BehlerFeatureFixture {};

  BOOST_AUTO_TEST_SUITE(behler_parinello_feature_tests);
  BOOST_AUTO_TEST_CASE(dummy) {}
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
