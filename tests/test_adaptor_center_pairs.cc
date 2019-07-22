/**
 * file   test_center_pairs.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jul 2018
 *
 * @brief  tests the implementation of the center pair adaptor
 *
 * Copyright  2018 Markus Stricker COSMO (EPFL), LAMMM (EPFL)
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

  // using Fixtures = boost::mpl::list<
  //     PairFixtureCenterPairs<PairFixtureSimple<StructureManagerCenters>>,
  //     PairFixtureCenterPairs<PairFixtureCenters>
  //     // PairFixtureCenterPairs<ManagerFixture<StructureManagerLammps>>
  //     >;

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
   * Iteration test for center pairs adaptor. It also checks if the first pair
   * in each atoms neighbourhood is itself.
   */
  //todo(markus) this filtering does not seem to work with ghost atoms
  BOOST_FIXTURE_TEST_CASE(iterator_test, PairFixtureCenters) {
                          //PairFixtureSimple<StructureManagerCenters>) {
    // auto adaptor_center_pairs{
    //     make_adapted_manager<AdaptorCenterPairs>(pair_manager, cutoff)};
    // adaptor_center_pairs->update();

    auto adaptor_center_pairs{
        make_adapted_manager<AdaptorCenterPairs>(pair_manager, cutoff)};
    adaptor_center_pairs->update();

    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{true};

    for (auto atom : pair_manager) {
      auto index{atom.get_global_index()};
      std::cout << "pairpair index " << index << std::endl;
      BOOST_CHECK_EQUAL(index, atom_counter);
      ++atom_counter;

      auto type{atom.get_atom_type()};

      for (auto pair : atom) {
        auto pair_offset{pair.get_global_index()};
        auto pair_type{pair.get_atom_type()};
        if (verbose) {
          std::cout << "pairpair (" << atom.back() << ", " << pair.back()
                    << "), pair_counter = " << pair_counter
                    << ", pair_offset = " << pair_offset
                    << ", atom types = " << type << ", " << pair_type
                    << std::endl;
        }

        BOOST_CHECK_EQUAL(pair_counter, pair_offset);
        ++pair_counter;
      }
    }

    int ctr{0};
    atom_counter = 0;
    pair_counter = 0;

    std::cout << "ACP-test test size() " << adaptor_center_pairs->size() << std::endl;
    std::cout << "ACP-test test get_size() " << adaptor_center_pairs->get_size()
              << std::endl;
    std::cout << "ACP-test test get_size_with_ghosts() "
              << adaptor_center_pairs->get_size_with_ghosts() << std::endl;

    for (auto atom : adaptor_center_pairs) {
      std::cout << "ctr " << ++ctr << std::endl;

      auto index{atom.get_global_index()};
      std::cout << "index atom " << index << std::endl;
      BOOST_CHECK_EQUAL(index, atom_counter);
      ++atom_counter;

      auto type{atom.get_atom_type()};

      for (auto pair : atom) {
        auto pair_offset{pair.get_global_index()};
        auto pair_type{pair.get_atom_type()};
        if (verbose) {
          std::cout << "index pair (" << atom.back() << ", " << pair.back()
                    << "), pair_counter = " << pair_counter
                    << ", pair_offset = " << pair_offset
                    << ", atom types = " << type << ", " << pair_type
                    << std::endl;
        }

        BOOST_CHECK_EQUAL(pair_counter, pair_offset);
        ++pair_counter;
      }
    }
    auto natoms_current{adaptor_center_pairs->get_size()};
    std::cout << "natoms_current " << natoms_current << std::endl;
    auto natoms_parent{this->pair_manager->get_size()};
    std::cout << "natoms_parent " << natoms_parent << std::endl;

    auto pairs_without_ii{this->pair_manager->get_nb_clusters(2)};
    std::cout << "pairs_without_ii " << pairs_without_ii << std::endl;
    auto pairs_with_ii{adaptor_center_pairs->get_nb_clusters(2)};
    std::cout << "pairs_with_ii " << pairs_with_ii << std::endl;
    std::cout << "atom_counter " << atom_counter << std::endl;


    std::cout << "Number of atoms parent " << natoms_parent << std::endl;
    std::cout << "Number of pairs parent " << pairs_without_ii << std::endl;
    std::cout << "Number of pairs current " << pairs_with_ii << std::endl;
    // Check increase in number of pairs (should be number of atoms, because
    // only ii-pairs are added
    BOOST_CHECK_EQUAL(pairs_without_ii + natoms_parent, pairs_with_ii);
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
