/**
 * file   test_center_pairs.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jul 2018
 *
 * @brief  tests the implementation of the center pair adaptor
 *
 * Copyright  2018 Till Junge, Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"

#include <vector>

namespace rascal {

  using Fixtures = boost::mpl::list<
      PairFixtureCenterPairs<PairFixtureSimple<StructureManagerCenters>>,
      PairFixtureCenterPairs<PairFixtureCenters>
      // PairFixtureCenterPairs<ManagerFixture<StructureManagerLammps>>
      >;

  BOOST_AUTO_TEST_SUITE(center_pairs_adaptor_test);

  /* ---------------------------------------------------------------------- */
  /**
   * test center pair constructor, this is done here instead
   * of in the Fixture, because of the default constructor is deleted.
   */
  BOOST_FIXTURE_TEST_CASE(constructor_test, PairFixtureCenters) {
    auto adaptor_center_pairs{
        make_adapted_manager<AdaptorCenterPairs>(pair_manager, cutoff)};
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Update test
   */
  BOOST_FIXTURE_TEST_CASE(update_test, PairFixtureCenters) {
    auto adaptor_center_pairs{
        make_adapted_manager<AdaptorCenterPairs>(pair_manager, cutoff)};
    adaptor_center_pairs->update();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Iteration test for strict adaptor. It also checks if the first pair in each
   * atoms neighbourhood is itself.
   */
  BOOST_FIXTURE_TEST_CASE(iterator_test, PairFixtureCenterPairs) {
    // auto adaptor_center_pairs{
    //     make_adapted_manager<AdaptorCenterPairs>(pair_manager, cutoff)};
    // adaptor_center_pairs->update();

    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{false};

    for (auto atom : this->adaptor_center_pairs) {
      auto index{atom.get_global_index()};
      BOOST_CHECK_EQUAL(index, atom_counter);

      auto type{atom.get_atom_type()};
      ++atom_counter;

      for (auto pair : atom) {
        auto pair_offset{pair.get_global_index()};
        auto pair_type{pair.get_atom_type()};
        if (verbose) {
          std::cout << "pair (" << atom.back() << ", " << pair.back()
                    << "), pair_counter = " << pair_counter
                    << ", pair_offset = " << pair_offset
                    << ", atom types = " << type << ", " << pair_type
                    << std::endl;
        }

        BOOST_CHECK_EQUAL(pair_counter, pair_offset);
        ++pair_counter;
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
