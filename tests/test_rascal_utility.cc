/**
 * file   test_rascal_utility.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   14 Oct 2019
 *
 * @brief Tests the implementation of interpolator
 *
 * Copyright © 2019  Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
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

#include "tests.hh"
#include "rascal_utility.hh"
#include "representations/calculator_spherical_expansion.hh"

namespace rascal {

  // TODO(all) not sure about the naming convention of tests, camelcase ...
  BOOST_AUTO_TEST_SUITE(RascalUtilityTests);

  using internal::RadialBasisType;
  using internal::AtomicSmearingType;
  using internal::OptimizationType;

  /*
   * enum class RadialBasisType { GTO, DVR, End_ };
   * enum class AtomicSmearingType { Constant, PerSpecies, Radial, End_ };
   * enum class OptimizationType { None, Interpolator, End_ };
   *
   * enum value tuple  ->  combined key
   *    (0,0,0)                 0
   *    (1,0,0)                 1
   *    (0,1,0)                 2
   *    (1,1,0)                 3
   *    (0,2,0)                 4
   *    (1,2,0)                 5
   *    (0,0,1)                 6
   *    (1,0,1)                 7
   *    (0,1,1)                 8
   *    (1,1,1)                 9
   *    (0,2,1)                10
   *    (1,2,1)                11
   */
  BOOST_AUTO_TEST_CASE(combine_enums_test) {
    BOOST_CHECK_EQUAL(combine_enums(RadialBasisType::GTO,
                                    AtomicSmearingType::Constant,
                                    OptimizationType::None),
                      0);

    BOOST_CHECK_EQUAL(combine_enums(RadialBasisType::DVR,
                                    AtomicSmearingType::PerSpecies,
                                    OptimizationType::None),
                      3);
    BOOST_CHECK_EQUAL(combine_enums(RadialBasisType::DVR,
                                    AtomicSmearingType::Constant,
                                    OptimizationType::Interpolator),
                      7);
    BOOST_CHECK_EQUAL(combine_enums(RadialBasisType::DVR,
                                    AtomicSmearingType::Radial,
                                    OptimizationType::Interpolator),
                      11);
  }
  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
