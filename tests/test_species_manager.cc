/**
 * file    test_species_manager.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   14 Sep 2018
 *
 * @brief tests the implementation of the adaptor filter species
 *
 * Copyright  2018 Markus Stricker, Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "structure_managers/species_manager.hh"
#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(adaptor_filter_species_test);

  template <class ManagerImplementation, size_t MaxOrder>
  struct SpeciesManagerFixture {
    using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;

    SpeciesManagerFixture() : species_manager{fixture.manager} {}

    ~SpeciesManagerFixture() = default;

    size_t get_MaxOrder() { return MaxOrder; }

    ManagerFixture<ManagerImplementation> fixture{};
    SpeciesManager_t species_manager;
  };

  using Fixtures =
      boost::mpl::list<SpeciesManagerFixture<StructureManagerCenters, 1>,
                       SpeciesManagerFixture<StructureManagerLammps, 2>>;

  using FixturesMax1 =
      boost::mpl::list<SpeciesManagerFixture<StructureManagerCenters, 1>>;
  using FixturesMax2 =
      boost::mpl::list<SpeciesManagerFixture<StructureManagerLammps, 2>>;
  using FixturesMax3 = boost::mpl::list<>;

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Fixtures, Fix) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_species_test, Fix, FixturesMax1, Fix) {
    std::map<int, int> species_counter{};
    for (auto && atom : Fix::fixture.manager) {
      species_counter[atom.get_atom_type()]++;
    }

    Fix::species_manager.update();

    for (auto && tup : species_counter) {
      auto species{tup.first};
      auto nb_atoms{tup.second};
      auto nb_filtered{
          Fix::species_manager[std::array<int, 1>{{species}}].size()};
      BOOST_CHECK_EQUAL(nb_atoms, nb_filtered);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(pair_species_test, Fix, FixturesMax2, Fix) {
    std::map<std::array<int, 2>, int> species_counter{};
    for (auto && atom : Fix::fixture.manager) {
      for (auto && pair : atom) {
        species_counter[pair.get_atom_types()]++;
      }
    }

    Fix::species_manager.update();
    for (auto && tup : species_counter) {
      auto species{tup.first};
      auto nb_clusters{tup.second};
      auto nb_filtered{Fix::species_manager[species].nb_clusters(2)};
      BOOST_CHECK_EQUAL(nb_clusters, nb_filtered);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
