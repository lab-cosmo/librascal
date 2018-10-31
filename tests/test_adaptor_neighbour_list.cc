/**
 * file   test_adaptor_neighbour_list.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Oct 2018
 *
 * @brief tests the implementation of the adaptor for building a
 * neighbour list, depends on traits if it is full of minimal
 *
 * Copyright © 2018 Markus Stricker, Till Junge, COSMO (EPFL), LAMMM (EPFL)
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
  /*
   * very simple 9 atom neighbour list build without periodicity
   *
   * ``manager`` is a MaxOrder=1 atom manager
   */
  BOOST_FIXTURE_TEST_CASE(simple_cubic_9_neighbour_list,
                          ManagerFixtureFile<StructureManagerCenters>) {

    constexpr bool verbose{false};

    AdaptorNeighbourList<StructureManagerCenters> SM2{manager, cutoff};
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

  /* ---------------------------------------------------------------------- */
  //! test if hcp managers are constructed
  BOOST_FIXTURE_TEST_CASE(constructor_test_hcp,
                          ManagerFixtureNeighbourTwoHcp) {
  }

  /* ---------------------------------------------------------------------- */
  //! test if fcc managers are constructed
  BOOST_FIXTURE_TEST_CASE(constructor_test_fcc,
                          ManagerFixtureNeighbourCheckFcc
                          <StructureManagerCenters>) {
  }

  /* ---------------------------------------------------------------------- */
  /**
   * simple neighbourhood test with periodicity only in x-direction and a check
   * for internal consistency
   *
   * ``manager`` is a StructureManager with MaxOrder=1
   */
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
  /*
   * test if two differently defined 2-atom units cells of hcp crystal structure
   * yield the same number of neighbours per atom, if the cutoff is increased.
   *
   * ``manager_1`` and ``manager_2`` each hold a different unit cell for a hcp
   * crystal system.
   */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test_hcp,
                          ManagerFixtureNeighbourTwoHcp) {

    /*
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
  /*
   * Test if two differently defined 1-atom and 4-atom units cells of fcc
   * crystal structure yield the same number of neighbours per atom at the
   * origin, if the cutoff is increased. ``manager_1`` has one atom at the
   * origin,``manager_2`` has 4 atoms (tetraeder). This is done to check if
   * skewed-ness affects the neighbour list algorithm.
   *
   * ``manager_1`` and ``manager_2`` each hold a different unit cell for a fcc
   * crystal system.
   */
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
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          bool is_in{dist < cutoff_tmp};
          if (verbose) {
            std::cout << "1 pair ("
                      << atom.get_position().transpose() << ", "
                      << pair.get_position().transpose() << ", "
                      << dist << ", " << cutoff_tmp << ", " << is_in
                      << ")" <<std::endl;
          }
          if (is_in) {
            neighbours_per_atom1.back()++;
          }
        }
      }

      for (auto atom : pair_manager2) {
        neighbours_per_atom2.push_back(0);
        for (auto pair : atom) {
          double dist = {(atom.get_position()
                          - pair.get_position()).norm()};
          bool is_in{dist < cutoff_tmp};
          if (verbose) {
            std::cout << "2 pair ("
                      << atom.get_position().transpose() << ", "
                      << pair.get_position().transpose() << ", "
                      << dist << ", " << cutoff_tmp << ", " << is_in
                      << ")" <<std::endl;
          }
          if (is_in) {
            neighbours_per_atom2.back()++;
          }
        }
      }

      /*
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

  /* ---------------------------------------------------------------------- */
  /*
   * Test if an increasingly skewed unit cell gives the same number of
   * neighbours than the non-skewed one. Skewing is achieved by the
   * multiplier.
   *
   * The fixture does not provide a manager. Instead, it provies basic cubic
   * unit cell, positions and the shear transformation ``skew_multiplier`` for
   * the unit cell. This is used in the loop to construct increasingly skewed
   * unit cells; positions are shifted accordingly to stay inside the new unit
   * cell.
   */
  BOOST_FIXTURE_TEST_CASE(test_neighbour_list_skewed,
                          ManagerFixtureSkew) {

    constexpr static bool verbose{false};

    int ncells{3};

    // container for storing the number of neighbours per atom in all cells
    std::vector<std::vector<int>> neighbours;
    neighbours.resize(ncells);

    // loop over cells
    for (int i{0}; i < ncells; ++i) {
      // manager constructed within this loop
      StructureManagerCenters manager;

      if (verbose) std::cout << "------------ cells " << i << std::endl;

      // get reference data
      auto skewer{skew_multiplier};
      // change shear multiplier
      skewer(0,1) *= i;
      // calculate unit cell
      auto cell_skw = skewer * cell;
      auto cell_skw_inv{cell_skw.inverse().eval()};

      if (verbose) {
        std::cout << "cell and inverse " << std::endl;
        std::cout << cell_skw << std::endl;
        std::cout << cell_skw_inv << std::endl;
      }

      // change initial atomic positions according to skewedness
      auto pos_skw{positions};
      for (int j{0}; j < natoms; ++j) {
        auto p = pos_skw.col(j);
        if (verbose) {
          std::cout << " >>>>>>> original position atom " << j << "\n"
                    << p.transpose() << std::endl;
        }
        // find minimal multiplier to project position out of cubic unit cell
        auto mult{(cell_skw_inv * p).eval()};
        auto mult_floor{mult};
        for (auto m{0}; m < 3; ++m) {mult_floor[m] = std::floor(mult[m]);}
        auto m{mult_floor.cwiseAbs().eval()};

        if (verbose) {
          std::cout << " mult \n" << mult.transpose() << std::endl;
          std::cout << "  mult_floor \n" << m.transpose() << std::endl;
          std::cout << "  skewed position  <<<<<<<\n" << p.transpose()
                    << std::endl;
        }
        // set new position
        pos_skw.col(j) = p + cell_skw * m;
      }

      if (verbose) {
        std::cout << "positions skewed \n" << pos_skw << std::endl;
      }

      // construct manager with skewed unit cell and shifted positions
      manager.update(pos_skw, atom_types, cell_skw,
                     Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});

      // build neighbourlist
      AdaptorNeighbourList<StructureManagerCenters> pair_manager{manager,
          cutoff};
      pair_manager.update();

      // make strict for counting neighbours
      AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>
        adaptor_strict{pair_manager, cutoff};
      adaptor_strict.update();

      // count strict neighbours
      for (auto atom : adaptor_strict) {
        neighbours[i].push_back(0);
        for (auto pair : atom) {
          neighbours[i].back()++;
        }
      }
    }
    for (auto i{0}; i < ncells; ++i) {

    }
    // check neighbours
    for (auto i{1}; i < ncells; ++i) {
      if (verbose) {
        std::cout << "cell " << i << std::endl;
        for (auto neigh_per_atom : neighbours[i]) {
          std::cout << neigh_per_atom << " ";
        }
        std::cout << std::endl;
      }

      // neighbours[i] are the skewed unit cells with adapted positions
      BOOST_CHECK_EQUAL_COLLECTIONS(neighbours[i].begin(),
                                    neighbours[i].end(),
                                    // expected values from unskewed cell
                                    neighbours[0].begin(),
                                    neighbours[0].end());
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
