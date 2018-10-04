/**
 * file   test_adaptor_half_neighbour_list.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief  tests the implementation of the half neighbourlist adaptor
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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
#include "test_structure.hh"
#include "structure_managers/adaptor_half_neighbour_list.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE(half_neighbourlist_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerLammps>) {

    int npairs_full{0};
    for (auto atom : manager) {
      for (auto pair : atom) {
        npairs_full++;
      }
    }


    std::cout << "Setting up strict manager" << std::endl;
    AdaptorHalfList<StructureManagerLammps> adaptor{manager};
    adaptor.update();

    int npairs_half{0};
    // std::cout << "Testing adaptor_strict" << std::endl;
    for (auto atom : adaptor) {
      // std::cout << "strict atom out " << atom.get_index() << std::endl;
      for (auto pair : atom) {
        npairs_half++;
      }
    }

    std::cout << "Full/half " << npairs_full << "/" << npairs_half << std::endl;
  }




  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
