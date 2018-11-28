/**
 * file   test_adaptor_strict.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief  tests the implementation of the strict structure adaptor
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_structure.hh"
#include "structure_managers/adaptor_strict.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE(strict_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerCell>) {
    double cut_off{0.9*cutoff_max};

    // std::cout << "Testing manager cell iteration" << std::endl;
    // for (auto atom : manager) {
    //   for (auto pair : atom) {
    //     std::cout << "out " << pair.get_index() << std::endl;
    //   }
    // }

    // TODO: Check if the neighbour list is actually a linked cell list, not
    // just all neighbours

    // std::cout << "Setting up strict manager" << std::endl;
    AdaptorStrict<StructureManagerCell> adaptor{manager, cut_off};
    adaptor.update();

    // std::cout << "Testing adaptor_strict" << std::endl;
    for (auto atom : adaptor) {
      // std::cout << "strict atom out " << atom.get_index() << std::endl;
      for (auto pair : atom) {
        // TODO: commented out after Félix's repair
        // std::cout << "  strict pair out " << pair.get_index() << std::endl;
        // auto dist = adaptor.get_distance(pair);
        // std::cout << "distance " << dist << std::endl;
      }
    }
  }


  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
