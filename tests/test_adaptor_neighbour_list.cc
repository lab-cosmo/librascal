/**
 * file   test_adaptor_neighbour_list.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   05 Oct 2018
 *
 * @brief tests the implementation of the adaptor for building a
 * neighbour list, depends on traits if it is full of minimal
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
#include "structure_managers/adaptor_neighbour_list.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(neighbour_list_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_fcc,
                          ManagerFixtureNeighbourCheckFcc
                          <StructureManagerCenters>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_hcp,
                          ManagerFixtureNeighbourCheckHcp
                          <StructureManagerCenters>) {
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(test_build_neighbour_list_from_atoms,
                          ManagerFixtureSimple<StructureManagerCenters>){

    constexpr bool verbose{false};

    if (verbose) std::cout << "===> zeroth order manager " << std::endl;
    //! testing iteration of zerot-th order manager
    for (auto atom : manager) {
      if (verbose) {
        std::cout << "atom " << atom.back() << std::endl;
      }
    }

    if (verbose) std::cout << "<== zeroth order manager " << std::endl;

    AdaptorNeighbourList<StructureManagerCenters> pair_manager{manager, cutoff};
    pair_manager.update();

    auto n_pairs{0};
    for (auto atom : pair_manager) {
      if (verbose) std::cout << "atom " << atom.back() << std::endl;
      for (auto pair : atom) {
        n_pairs++;
        if (verbose) {
          std::cout << "   complete pair "
                    << atom.back() << " " << pair.back()
                    << " glob " << pair.get_global_index() << std::endl;
        }
      }
    }
    if (verbose) std::cout << "Number of pairs " << n_pairs << std::endl;
    BOOST_CHECK_EQUAL(n_pairs, pair_manager.get_nb_clusters(2));
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test_hcp,
                          ManagerFixtureNeighbourCheckHcp
                          <StructureManagerCenters>) {

    /**
     * Note: since the cell vectors are different, it is possible that one of
     * the two atoms is repeated into a different cell due to periodicity. This
     * leads to a difference in number of neighbours. Therefore the strict
     * cutoff is check to ensure the exakt same number of neighbours.
     */

    constexpr bool verbose{false};

    if(verbose) std::cout << "HCP test " << cutoff << std::endl;

    int mult = 10;

    for (auto i{1}; i < mult; ++i) {
      auto cutoff_tmp = i * cutoff;

      std::vector<int> neighbours_per_atom1{};
      std::vector<int> neighbours_per_atom2{};

      neighbours_per_atom1.resize(0);
      neighbours_per_atom1.resize(0);

      if (verbose) {
        std::cout << "hcp test cutoff " << cutoff_tmp << std::endl;
      }

      AdaptorNeighbourList<StructureManagerCenters> pair_manager1{manager_1,
          cutoff_tmp};
      pair_manager1.update();

      AdaptorNeighbourList<StructureManagerCenters> pair_manager2{manager_2,
          cutoff_tmp};
      pair_manager2.update();

      if (verbose) std::cout << "Manager 1" << std::endl;

      for (auto atom : pair_manager1) {
        neighbours_per_atom1.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "1 pair "
                      << atom.back() << " "
                      << pair.back() << std::endl;
          }
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom1.back()++;
          }
        }
      }

      if (verbose) std::cout << "Manager 2" << std::endl;

      for (auto atom : pair_manager2) {
        neighbours_per_atom2.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "2 pair "
                      << atom.back() << " "
                      << pair.back() << std::endl;
          }
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom2.back()++;
          }
        }
      }

      BOOST_CHECK_EQUAL_COLLECTIONS(neighbours_per_atom1.begin(),
                                    neighbours_per_atom1.end(),
                                    neighbours_per_atom2.begin(),
                                    neighbours_per_atom2.end());

      for (auto i{0}; i < natoms; ++i) {
        if (verbose) {
          std::cout << "neigh1/neigh2: i " << i << " "
                    << neighbours_per_atom1[i] << "/"
                    << neighbours_per_atom2[i] << std::endl;
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test_fcc,
                          ManagerFixtureNeighbourCheckFcc
                          <StructureManagerCenters>) {

    constexpr bool verbose{false};

    if (verbose) std::cout << "FCC test " << std::endl;

    int mult = 8;

    for (auto i{1}; i < mult; ++i) {
      auto cutoff_tmp = i * cutoff;

      std::vector<int> neighbours_per_atom1{};
      std::vector<int> neighbours_per_atom2{};

      neighbours_per_atom1.resize(0);
      neighbours_per_atom1.resize(0);

      if (verbose) {
        std::cout << "fcc cutoff " << cutoff_tmp << std::endl;
      }

      AdaptorNeighbourList<StructureManagerCenters> pair_manager1{manager_1,
          cutoff_tmp};
      pair_manager1.update();

      AdaptorNeighbourList<StructureManagerCenters> pair_manager2{manager_2,
          cutoff_tmp};
      pair_manager2.update();

      for (auto atom : pair_manager1) {
        neighbours_per_atom1.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "1 pair "
                      << atom.back() << " "
                      << pair.back() << std::endl;
          }
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom1.back()++;
          }
        }
      }

      for (auto atom : pair_manager2) {
        neighbours_per_atom2.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "2 pair "
                      << atom.back() << " "
                      << pair.back() << std::endl;
          }
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom2.back()++;
          }
        }
      }

      /**
       * only the first index atom can be checked, since the cell with only one
       * atom does not allow for comparison with other atom's number of
       * neighbours
       */
      BOOST_CHECK_EQUAL(neighbours_per_atom1[0],
                        neighbours_per_atom2[0]);
      if (verbose) {
        std::cout << "neigh1/neigh2: "
                  << neighbours_per_atom1[0] << "/"
                  << neighbours_per_atom2[0] << std::endl;
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
