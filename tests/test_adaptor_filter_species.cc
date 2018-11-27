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
#include "structure_managers/species_manager.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE(adaptor_filter_species_test);

  template<class ManagerImplementation, size_t MaxOrder>
  struct SpeciesManagerFixture
  {
    using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;

    SpeciesManagerFixture():
      species_manager{fixture.manager}{}

    ~SpeciesManagerFixture() = default;

    size_t get_MaxOrder() {return MaxOrder;}

    ManagerFixture<ManagerImplementation> fixture{};
    SpeciesManager_t species_manager;
  };

  using Fixtures = boost::mpl::list<
    SpeciesManagerFixture<StructureManagerCenters, 1>,
    SpeciesManagerFixture<StructureManagerLammps, 2> >;

  using FixturesMax1 = boost::mpl::list<
    SpeciesManagerFixture<StructureManagerCenters, 1>>;
  using FixturesMax2 = boost::mpl::list<
    SpeciesManagerFixture<StructureManagerLammps, 2>>;
  using FixturesMax3 = boost::mpl::list<>;

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Fixtures, Fix) {

//    constexpr bool verbose{true};
//
//    /* skeleton test */
//    std::cout << "Skeleton test" << std::endl;
//
//    for (auto atom : manager) {
//      if (verbose) {
//        std::cout << "Atom "
//                  << atom.back() << " type "
//                  << atom.get_atom_type()
//        	  << std::endl;
//      }
//    }
  }

  BOOST_AUTO_TEST_SUITE_END();



}  // rascal
