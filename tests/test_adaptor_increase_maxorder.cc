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
#include "test_structure.hh"
#include "structure_managers/adaptor_increase_maxorder.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(maxlevel_increase_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_order_zero,
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

    AdaptorMaxOrder<StructureManagerCenters> pair_manager{manager, cutoff};
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
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerChain>){

    constexpr bool verbose{false};
    constexpr bool check_below{false};

    // Check underlying manager
    if (check_below) std::cout << ">============ below" << std::endl;
    size_t npairs1{0};
    if (check_below) {
      for (auto atom : manager_chain) {
        if (verbose) {
          std::cout << "chain atom "
                    << atom.back()
                    << std::endl;
        }
        for (auto pair : atom) {
          npairs1++;
          if (verbose) {
            std::cout << " chain pair "
                      << pair.back()
                      << " glob " << pair.get_global_index()
                      << std::endl;
          }
        }
      }
    }
    if (check_below) {
      std::cout << "number of pairs " << npairs1 << std::endl;
      std::cout << "<============ below" << std::endl;
    }

    AdaptorMaxOrder<StructureManagerChain> adaptor{manager_chain, cutoff};
    adaptor.update();

    if (verbose) {
      std::cout << "Adaptor increase MaxOrder" << std::endl;
    }

    auto natoms{0};
    auto npairs{0};
    auto n_triplets{0};
    for (auto atom : adaptor) {
      natoms++;
      if (verbose) {
        std::cout << atom.back()
        	  << std::endl;
      }

      if (verbose) std::cout << "position: " << atom.get_position() << std::endl;

      for (auto pair : atom) {
        npairs++;
        if (verbose) {
          std::cout << "   complete pair "
                    << atom.back() << " " << pair.back()
                    << " glob " << pair.get_global_index() << std::endl;
        }
        for (auto triplet : pair) {
          n_triplets++;
          if (verbose) {
            std::cout << "             triplet "
                      << triplet.back()
                      << " global " << triplet.get_global_index()
                      << std::endl;
            std::cout << "                         complete "
                      << atom.back() << " "
                      << pair.back() << " "
                      << triplet.back() << std::endl;
          }
        }
      }
    }
    if(verbose) std::cout << "Number of triplets: " << n_triplets << std::endl;
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test,
                          ManagerFixtureNeighbourComparison
                          <StructureManagerCenters>) {

    constexpr bool verbose{false};

    std::vector<int> neighbours_per_atom1{};
    std::vector<int> neighbours_per_atom2{};

    neighbours_per_atom1.resize(0);
    neighbours_per_atom1.resize(0);

    AdaptorMaxOrder<StructureManagerCenters> pair_manager1{manager_1,
        cutoff};
    pair_manager1.update();

    AdaptorMaxOrder<StructureManagerCenters> pair_manager2{manager_2,
        cutoff};
    pair_manager2.update();

    for (auto atom : pair_manager1) {
      neighbours_per_atom1.push_back(0);
      for (auto pair : atom) {
        if (verbose) {
          std::cout << "1 pair "
                    << atom.back() << " "
                    << pair.back() << std::endl;
        }
        neighbours_per_atom1.back()++;
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
        neighbours_per_atom2.back()++;
      }
    }
    for (auto i{0}; i < natoms; ++i) {
      if (verbose) {
      std::cout << "neigh1/neigh2: "
                << neighbours_per_atom1[i] << "/"
                << neighbours_per_atom2[i] << std::endl;
      }
      BOOST_CHECK_EQUAL(neighbours_per_atom1[i],
                        neighbours_per_atom2[i]);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
