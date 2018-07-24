/**
 * file   test_adaptor_increase_maxlevel.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   20 Jun 2018
 *
 * @brief tests the implementation of the adaptor increase maxlevel
 * (atom list to pairs, pairs to triplets, etc.)
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
#include "test_neighbourhood.hh"
#include "neighbourhood_managers/adaptor_increase_maxlevel.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(maxlevel_increase_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test, ManagerFixture_chain){

    constexpr bool verbose{true};
    constexpr bool check_below{true};

    // Check underlying manager
    if (check_below) {
      for (auto atom : manager_chain) {
        if (verbose) {
          std::cout << "chain atom "
                    << atom.back()
                    << std::endl;
        }
        for (auto pair : atom) {
          if (verbose) {
            std::cout << " chain pair "
                      << pair.back()
                      << std::endl;
          }
        }
      }
    }

    AdaptorMaxLevel<NeighbourhoodManagerChain> adaptor{manager_chain, cutoff};
    adaptor.update();

    if (verbose) {
      std::cout << "Adaptor increase MaxLevel" << std::endl;
    }

    for (auto atom : adaptor) {
      if (verbose) {
        std::cout << "increase MaxLevel atom "
        	  << atom.back()
        	  << std::endl;
      }
      for (auto pair : atom) {
        if (verbose) {
          std::cout << " increase MaxLevel pair "
        	    << pair.back()
        	    << std::endl;
        }
        std::cout << " increase MaxLevel before triplet " << std::endl;
        for (auto triplet : pair) {
          if (verbose) {
            std::cout << "    triplet "
                      << triplet.back()
                      << std::endl;
          }
        }
      }
    }
  }


  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
