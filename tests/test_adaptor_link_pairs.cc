/**
 * @file   test_link_pairs.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   28 Jan 2020
 *
 * @brief  tests the implementation of the swaping of pairs
 *
 * Copyright  2020 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#include "test_adaptor.hh"
#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>

namespace rascal {

  using multiple_fixtures = boost::mpl::list<
      MultipleStructureFixture<MultipleStructureManagerNLCCStrictFixture>>;

  BOOST_AUTO_TEST_SUITE(test_link_pairs);

  /* ---------------------------------------------------------------------- */
  /**
   * test strict neighbourhood constructor, this is done here instead
   * of in the Fixture, because of the default constructor is deleted.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, multiple_fixtures,
                                   Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      auto adaptor{make_adapted_manager<AdaptorLinkPairs>(manager)};
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test update
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_test, Fix, multiple_fixtures, Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      auto adaptor{make_adapted_manager<AdaptorLinkPairs>(manager)};
      adaptor->update();
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Iteration test for strict link pairs. It also checks if the types of the
   * original atoms are accessed correctly.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(iterator_test, Fix, multiple_fixtures, Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      auto adaptor{make_adapted_manager<AdaptorLinkPairs>(manager)};
      adaptor->update();

      auto structure = extract_underlying_manager<0>(adaptor);
      auto atom_types = structure->get_atom_types();

      int atom_counter{};
      int pair_counter{};
      constexpr bool verbose{false};

      for (auto atom : adaptor) {
        auto index{atom.get_global_index()};
        BOOST_CHECK_EQUAL(index, atom_counter);

        auto type{atom.get_atom_type()};
        BOOST_CHECK_EQUAL(type, atom_types[index]);
        ++atom_counter;

        for (auto pair : atom.pairs_with_self_pair()) {
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
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test that the atom index from a neighbour matches the atom tag of the
   * ClusterRefKey returned by get_atom_j for a manager using as root
   * implementation `StructureManagerCenters`.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(get_atom_j_test, Fix, multiple_fixtures,
                                   Fix) {
    auto && managers = Fix::managers;
    constexpr bool verbose{false};
    for (auto & manager : managers) {
      auto adaptor{make_adapted_manager<AdaptorLinkPairs>(manager)};
      adaptor->update();

      if (verbose) {
        std::cout << "adaptor size : " << adaptor->size() << std::endl;
        std::cout << "adaptor size_wg : " << adaptor->get_size_with_ghosts()
                  << std::endl;
      }
      for (auto atom : adaptor) {
        for (auto pair : atom.pairs()) {
          auto atom_j_index = adaptor->get_atom_index(pair.back());
          auto atom_j = pair.get_atom_j();
          auto atom_j_tag = atom_j.get_atom_tag_list();
          if (verbose) {
            std::cout << "neigh: " << atom_j_index
                      << " tag_j: " << atom_j_tag[0] << std::endl;
          }

          BOOST_CHECK_EQUAL(atom_j_index, atom_j_tag[0]);
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test the main functionality of AdaptorLinkPairs:
   *   + d_ij are equal d_ji and v_ij = - v_ji
   *   + proper access to v_ji and d_ji is maintained when other adaptors are
   *     stacked on top of AdaptorLinkPairs
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_swap_distance_test, Fix,
                                   multiple_fixtures, Fix) {
    bool verbose{false};
    auto & managers = Fix::managers;

    for (auto & pair_manager : managers) {
      auto adaptor{make_adapted_manager<AdaptorLinkPairs>(pair_manager)};
      adaptor->update();

      for (auto center : adaptor) {
        for (auto pair_ij : center.pairs_with_self_pair()) {
          auto && d_ij = adaptor->get_distance(pair_ij);
          auto && v_ij = adaptor->get_direction_vector(pair_ij);

          auto pair_ji = pair_ij.get_pair_ji();
          auto && d_ji = adaptor->get_distance(pair_ji);
          auto && v_ji = adaptor->get_direction_vector(pair_ji);
          BOOST_TEST(d_ij - d_ji < 1e-13);
          BOOST_TEST(((v_ij + v_ji).array() < 1e-13).all());
          if (verbose) {
            std::cout << "d_ij=" << d_ij << " d_ji=" << d_ji;
            std::cout << "v_ij=" << v_ij.transpose()
                      << " v_ji=" << v_ji.transpose() << std::endl;
          }
        }
      }

      // test access with aditional adaptors
      auto adaptor_triplets{make_adapted_manager<AdaptorMaxOrder>(adaptor)};
      adaptor_triplets->update();

      for (auto center : adaptor_triplets) {
        for (auto pair_ij : center.pairs_with_self_pair()) {
          auto && d_ij = adaptor_triplets->get_distance(pair_ij);
          auto && v_ij = adaptor_triplets->get_direction_vector(pair_ij);

          auto pair_ji = pair_ij.get_pair_ji();
          auto && d_ji = adaptor_triplets->get_distance(pair_ji);
          auto && v_ji = adaptor_triplets->get_direction_vector(pair_ji);
          BOOST_TEST(d_ij - d_ji < 1e-13);
          BOOST_TEST(((v_ij + v_ji).array() < 1e-13).all());
          if (verbose) {
            std::cout << "d_ij=" << d_ij << " d_ji=" << d_ji;
            std::cout << "v_ij=" << v_ij.transpose()
                      << " v_ji=" << v_ji.transpose() << std::endl;
          }
        }
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
