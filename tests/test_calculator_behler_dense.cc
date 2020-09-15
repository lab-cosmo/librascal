/**
 * @file   test_calculator_behler_dense.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@ruhr-uni-bochum.de>
 *
 * @date   15 Sep 2020
 *
 * @brief  tests for the mono-species-combination Behler-Parrinello calculator
 *
 * Copyright © 2020 Till Junge, Markus Stricker
 *
 * libRascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * libRascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libRascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it
 * with proprietary FFT implementations or numerical libraries, containing parts
 * covered by the terms of those libraries' licenses, the licensors of this
 * Program grant you additional permission to convey the resulting work.
 *
 */

#include "rascal/representations/calculator_behler_parrinello_dense.hh"
#include "rascal/utils/json_io.hh"

#include "behler_fixtures.hh"
#include "test_structure.hh"

#include <boost/test/unit_test.hpp>

#include <vector>

namespace rascal {

  struct CalculatorBehlerParinelloDenseFixture {
    CalculatorBehlerParinelloDenseFixture() : raw_params{json_io::load(filename)} {}
    std::string filename {};
    json raw_params;
  }

}  // rascal
