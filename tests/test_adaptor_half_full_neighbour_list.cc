/**
 * @file   test_adaptor_half_neighbour_list.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief tests the implementation of the half and full neighbourlist adaptors
 *
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/adaptor_full_neighbour_list.hh"
#include "structure_managers/adaptor_half_neighbour_list.hh"
#include "test_structure.hh"

#include <boost/test/unit_test.hpp>

constexpr double TOLERANCE = 1e-14;
namespace rascal {

  BOOST_AUTO_TEST_SUITE(half_neighbourlist_adaptor_test);
  /* ---------------------------------------------------------------------- */
  /*
   * test the reduction of a full to a half neighbour list
   *
   * ``manager`` is an object of type StructureManagerLammps, it has 3 atoms and
   * a full neighbour list (MaxOrder=2).
   */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerLammps>) {
    auto adaptor = make_adapted_manager<AdaptorHalfList>(manager);
    adaptor->update();
  }

  /* ---------------------------------------------------------------------- */
  /*
   * test if sum of distances are the same in the full list and 2 x the half
   * list
   *
   * ``manager`` is an object of type StructureManagerLammps, which has
   * MaxOrder=2 and a full neighbourlist.
   */
  BOOST_FIXTURE_TEST_CASE(iteration_and_distance_half_list,
                          ManagerFixture<StructureManagerLammps>) {
    constexpr bool verbose{false};

    double distance_sum_full{0.};

    int npairs_full{0};
    for (auto atom : manager) {
      for (auto pair : atom) {
        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_full += dist;
        npairs_full++;
      }
    }

    if (verbose) {
      std::cout << "Setting up half neighbourlist manager" << std::endl;
    }

    auto adaptor = make_adapted_manager<AdaptorHalfList>(manager);
    adaptor->update();

    double distance_sum_half{0.};

    int npairs_half{0};
    for (auto atom : adaptor) {
      if (verbose) {
        std::cout << "type " << atom.get_atom_type() << std::endl;
      }
      for (auto pair : atom) {
        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_half += dist;

        dist = {(pair.get_position() - atom.get_position()).norm()};
        distance_sum_half += dist;
        npairs_half++;
      }
    }

    if (verbose) {
      std::cout << "Full/half " << npairs_full << "/" << npairs_half
                << std::endl;
    }

    auto val{distance_sum_full - distance_sum_half};
    auto relative_error = val * val / (distance_sum_full * distance_sum_full);

    BOOST_CHECK_EQUAL(npairs_full, 4);
    BOOST_CHECK_EQUAL(npairs_half, 2);
    BOOST_CHECK(relative_error < TOLERANCE * TOLERANCE);
  }

  /* ---------------------------------------------------------------------- */
  /*
   * This test is checks if the reduction of a full to a half neighbourlist and
   * successive extension to a full one again leads to the same
   * results. ``manager`` is provided as an input from the fixture as a basis
   * for the adaptions.
   */
  BOOST_FIXTURE_TEST_CASE(full_to_half_to_full_neighbour_list,
                          ManagerFixture<StructureManagerLammps>) {
    constexpr bool verbose{false};

    // variable for summation for test position differences
    double distance_sum_full{0.};
    double distance_sum_full_half_full{0.};

    // half list construction
    auto adaptor_half = make_adapted_manager<AdaptorHalfList>(manager);
    adaptor_half->update();

    // back to full list again
    auto adaptor_full = make_adapted_manager<AdaptorFullList>(adaptor_half);
    adaptor_full->update();

    // iterate over initial manager with full list
    if (verbose) {
      std::cout << "---full---" << std::endl;
    }
    int npairs{0};
    for (auto atom : manager) {
      for (auto pair : atom) {
        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_full += dist;
        npairs++;
        if (verbose) {
          auto pair_offset{pair.get_global_index()};
          std::cout << "pair (" << atom.back() << ", " << pair.back() << "), "
                    << pair_offset << std::endl;
        }
      }
    }

    // iterate over half list adator (only for output, if verbose to compare
    // pair indices)
    if (verbose) {
      std::cout << "---full/half---" << std::endl;
    }
    for (auto atom : adaptor_half) {
      for (auto pair : atom) {
        if (verbose) {
          auto pair_offset{pair.get_global_index()};
          std::cout << "pair (" << atom.back() << ", " << pair.back() << "), "
                    << pair_offset << std::endl;
        }
      }
    }

    // iterate over extended half list adaptor
    if (verbose) {
      std::cout << "---full/half/full---" << std::endl;
    }
    int npairs_adapted{0};
    for (auto atom : adaptor_full) {
      for (auto pair : atom) {
        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_full_half_full += dist;
        npairs_adapted++;
        if (verbose) {
          auto pair_offset{pair.get_global_index()};
          std::cout << "pair (" << atom.back() << ", " << pair.back() << "), "
                    << pair_offset << std::endl;
        }
      }
    }

    auto val{distance_sum_full - distance_sum_full_half_full};
    auto relative_error = val * val / (distance_sum_full * distance_sum_full);
    // check if the sum and square of all distances is the same
    BOOST_CHECK(relative_error < TOLERANCE * TOLERANCE);

    // check counted number of pairs during iteration
    BOOST_CHECK_EQUAL(npairs, npairs_adapted);
    // check for same number of atoms
    BOOST_CHECK_EQUAL(manager->size(), adaptor_full->size());
    // check for number of atoms
    BOOST_CHECK_EQUAL(manager->get_nb_clusters(1),
                      adaptor_full->get_nb_clusters(1));
    // check for number of pairs
    BOOST_CHECK_EQUAL(manager->get_nb_clusters(2),
                      adaptor_full->get_nb_clusters(2));
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
