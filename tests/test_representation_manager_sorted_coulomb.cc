/**
 * file   test_representation_manager.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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
#include "test_representation_manager.hh"


namespace rascal {
  BOOST_AUTO_TEST_SUITE(representation_sorted_coulomb_test);
  /* ---------------------------------------------------------------------- */
  // BOOST_FIXTURE_TEST_CASE(constructor_test,
  // RepresentationFixture<StructureManagerJson>)
  // { 
    
  //   AdaptorNeighbourList<StructureManagerJson> nl{manager_json,cutoff};
  //   nl.update();
  //   AdaptorStrict<AdaptorNeighbourList<
  //                             StructureManagerJson>> strict_nl{nl,cutoff};
  //   strict_nl.update();

  //   using Representation_t = RepresentationManagerSortedCoulomb<
  //                  AdaptorStrict<AdaptorNeighbourList<StructureManagerJson>>>;

  //   Representation_t representation(strict_nl,central_decay,
  //                                   interaction_cutoff,interaction_decay,size);
  //   representation.compute();


  // }


  BOOST_AUTO_TEST_SUITE_END();
} // RASCAL