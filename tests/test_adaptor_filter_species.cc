/**
 * file    test_adaptor_filter_species.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   14 Sep 2018
 *
 * @brief tests the implementation of the adaptor filter species
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
#include "structure_managers/adaptor_filter_species.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(adaptor_filter_species_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_order_zero,
                          ManagerFixtureSimple<StructureManagerCenters>){

    constexpr bool verbose{true};

    /* skeleton test */
    std::cout << "Skeleton test" << std::endl;

    for (auto atom : manager) {
      if (verbose) {
        std::cout << "Atom "
                  << atom.back() << " type "
                  << atom.get_atom_type()
        	  << std::endl;
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();



}  // rascal
