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
#include "structure_managers/adaptor_half_neighbour_list.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(maxlevel_increase_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerChain>) {
    AdaptorMaxOrder<StructureManagerChain> adaptor{manager_chain};
    adaptor.update();
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iterator_test,
                          ManagerFixture<StructureManagerChain>) {

    constexpr bool verbose{false};
    constexpr bool check_below{false};

    // Check underlying manager
    if (check_below) std::cout << ">> underlying manager " << std::endl;
    size_t npairs1{0};
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
    if (check_below) {
      std::cout << "number of pairs " << npairs1 << std::endl;
      std::cout << "<< underlying manager" << std::endl;
    }

    auto npairs_tmp = manager_chain.get_nb_clusters(2);
    BOOST_CHECK_EQUAL(npairs_tmp, npairs1);

    AdaptorMaxOrder<StructureManagerChain> adaptor{manager_chain};
    adaptor.update();

    //! make sure the number of pairs gets carried over to the next layer
    auto npairs_adaptor = adaptor.get_nb_clusters(2);
    BOOST_CHECK_EQUAL(npairs_adaptor, npairs_tmp);

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

      if (verbose) std::cout << "position: " << atom.get_position
                     () << std::endl;

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

    // TODO: check for consistency in number of tuples
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test with 3 atoms, included stacking: full pair list -> half pair list ->
   * triplet list; SM is used as a shorthand for StructureManager. Checked
   * positions are specific to StructureManagerLammps.
   */
  BOOST_FIXTURE_TEST_CASE(pair_to_triplet_extension,
                          ManagerFixture<StructureManagerLammps>) {

    constexpr bool verbose{false};

    if (verbose) std::cout << ">> pair to triplet extension" << std::endl;

    // TODO: should this be in a fixture?
    AdaptorHalfList<StructureManagerLammps> SM2{manager};
    SM2.update();
    AdaptorMaxOrder<AdaptorHalfList<StructureManagerLammps>> SM3{SM2};
    SM3.update();

    // make sure number of pairs are carried over, since they are are not changed
    BOOST_CHECK_EQUAL(SM2.get_nb_clusters(2), SM3.get_nb_clusters(2));

    // only one possible triplet in this case?
    BOOST_CHECK_EQUAL(SM3.get_nb_clusters(3), 1);

    for (auto atom : SM3) {
      auto atom_index = atom.get_atom_index();
      auto atom_type = atom.get_atom_type();
      BOOST_CHECK_EQUAL(atom_type, type[atom_index]);

      auto atom_position = atom.get_position();
      for (auto pair : atom) {
        auto pair_index = pair.get_atom_index();
        auto pair_type = pair.get_atom_type();
        BOOST_CHECK_EQUAL(pair_type, type[pair_index]);

        auto pair_position = pair.get_position();
        auto diff_pos_pair = (pair_position - atom_position).norm();
        BOOST_CHECK_CLOSE(diff_pos_pair, 1., tol);

        for (auto triplet : pair) {
          if (verbose) {
            std::cout << "triplet " << atom.back() << " "
                    << pair.back() << " " << triplet.back() << std::endl;
          }
          auto triplet_index = triplet.get_atom_index();
          auto triplet_type = triplet.get_atom_type();
          BOOST_CHECK_EQUAL(triplet_type, type[triplet_index]);

          auto triplet_position = triplet.get_position();
          auto diff_pos_triplet = (triplet_position - atom_position).norm();
          BOOST_CHECK_CLOSE(diff_pos_triplet, 1., tol);
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(simple_cubic_8_extension,
                          ManagerFixture<StructureManagerJson>) {

    constexpr bool verbose{false};

    AdaptorNeighbourList<StructureManagerJson> SM2{manager_json, cutoff};
    SM2.update();

    auto npairs = SM2.get_nb_clusters(2);

    if (verbose) std::cout << "npairs " << npairs << std::endl;

    int np{0};
    for (auto atom : SM2) {
      for (auto pair : atom) {
        np++;
      }
    }
    if (verbose) std::cout << "np " << np << std::endl;
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
