/**
 * file   test_adaptor_strict.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief  tests the implementation of the strict neighbourhood adaptor
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
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
#include "test_neighbourhood.hh"
#include "neighbourhood_managers/adaptor_strict.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE(strict_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<NeighbourhoodManagerCell>) {
    double cut_off{cutoff_max/2};
    AdaptorStrict<NeighbourhoodManagerCell> adaptor{manager, cut_off};
    adaptor.update();
  }


  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
