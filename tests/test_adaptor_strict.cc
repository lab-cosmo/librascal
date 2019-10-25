/**
 * @file   test_adaptor_strict.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief  tests the implementation of the strict structure adaptor
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

#include "test_adaptor.hh"
#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <vector>

namespace rascal {

  using Fixtures = boost::mpl::list<
      PairFixtureStrict<PairFixtureSimple<StructureManagerCenters>>,
      PairFixtureStrict<PairFixtureCenters>
      // PairFixtureStrict<ManagerFixture<StructureManagerLammps>>
      >;

  using multiple_fixtures = boost::mpl::list<
      MultipleStructureFixture<MultipleStructureManagerNLFixture>,
      MultipleStructureFixture<MultipleStructureManagerNLCCFixture>>;

  BOOST_AUTO_TEST_SUITE(adaptor_strict_test);

  /* ---------------------------------------------------------------------- */
  /**
   * test strict neighbourhood constructor, this is done here instead
   * of in the Fixture, because of the default constructor is deleted.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, multiple_fixtures,
                                   Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      double cutoff{manager->get_cutoff()};
      auto adaptor{make_adapted_manager<AdaptorStrict>(manager, cutoff)};
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test update
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(update_test, Fix, multiple_fixtures, Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      double cutoff{manager->get_cutoff()};
      auto adaptor{make_adapted_manager<AdaptorStrict>(manager, cutoff)};
      adaptor->update();
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Iteration test for strict adaptor. It also checks if the types of the
   * original atoms are accessed correctly.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(iterator_test, Fix, multiple_fixtures, Fix) {
    auto && managers = Fix::managers;
    for (auto & manager : managers) {
      double cutoff{manager->get_cutoff()};
      auto adaptor_strict{make_adapted_manager<AdaptorStrict>(manager, cutoff)};
      adaptor_strict->update();

      auto structure = extract_underlying_manager<0>(adaptor_strict);
      auto atom_types = structure->get_atom_types();

      int atom_counter{};
      int pair_counter{};
      constexpr bool verbose{false};

      for (auto atom : adaptor_strict) {
        auto index{atom.get_global_index()};
        BOOST_CHECK_EQUAL(index, atom_counter);

        auto type{atom.get_atom_type()};
        BOOST_CHECK_EQUAL(type, atom_types[index]);
        ++atom_counter;

        for (auto pair : atom.with_self_pair()) {
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
      double cutoff{manager->get_cutoff()};
      auto adaptor_strict{make_adapted_manager<AdaptorStrict>(manager, cutoff)};
      adaptor_strict->update();

      for (auto atom : adaptor_strict) {
        for (auto pair : atom) {
          auto atom_j_index = adaptor_strict->get_atom_index(pair.back());
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

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_strict_test, Fix, multiple_fixtures,
                                   Fix) {
    bool verbose{false};
    auto & managers = Fix::managers;

    for (auto & pair_manager : managers) {
      double cutoff{pair_manager->get_cutoff()};
      std::vector<std::vector<int>> neigh_ids{};
      std::vector<std::vector<double>> neigh_dist{};
      std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec{};
      std::vector<std::vector<int>> neigh_ids_strict{};
      std::vector<std::vector<double>> neigh_dist_strict{};
      std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec_strict{};

      if (verbose) {
        std::cout << "Setting up strict manager with rc = " << cutoff
                  << std::endl;
      }
      auto adaptor_strict{
          make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
      adaptor_strict->update();

      if (verbose) {
        std::cout << "Setting up comparison list with rc = " << cutoff
                  << std::endl;
      }
      for (auto center : pair_manager) {
        std::vector<int> indices{};
        std::vector<double> distances{};
        std::vector<std::array<double, 3>> dir_vecs{};

        if (verbose) {
          // get_index returns iteration index
          std::cout << "cell atom out " << center.get_index();
          // get_atom_tag returns index from
          std::cout << " " << center.get_atom_tag() << " ";

          for (int ii{0}; ii < 3; ++ii) {
            std::cout << center.get_position()[ii] << " ";
          }
          std::cout << " " << center.get_atom_type() << std::endl;
        }

        for (auto neigh : center) {
          double distance{
              (center.get_position() - neigh.get_position()).norm()};
          if (distance <= cutoff) {
            indices.push_back(neigh.get_atom_tag());
            distances.push_back(distance);
            auto dir_vec{
                (neigh.get_position() - center.get_position()).array() /
                distance};
            std::array<double, 3> aa{{dir_vec(0), dir_vec(1), dir_vec(2)}};
            dir_vecs.push_back(aa);
            if (verbose) {
              std::cout << "cell neigh out " << neigh.get_index();
              std::cout << " " << neigh.get_atom_tag() << " ";

              for (int ii{0}; ii < 3; ++ii) {
                std::cout << neigh.get_position()[ii] << " ";
              }
              std::cout << " " << neigh.get_atom_type() << std::endl;
            }
          }
        }
        neigh_ids.push_back(indices);
        neigh_dist.push_back(distances);
        neigh_dir_vec.push_back(dir_vecs);
        // break;
      }

      if (verbose) {
        std::cout << "Setting get adaptor_strict info" << std::endl;
      }
      for (auto center : adaptor_strict) {
        // auto icenter{center.get_index()};
        std::vector<int> indices_{};
        std::vector<double> distances_{};
        std::vector<std::array<double, 3>> dir_vecs_{};

        if (verbose) {
          // get_index returns iteration index
          std::cout << "strict atom out " << center.get_index();
          // get_atom_tag returns index from
          std::cout << " " << center.get_atom_tag() << " ";

          for (int ii{0}; ii < 3; ++ii) {
            std::cout << center.get_position()[ii] << " ";
          }
          std::cout << " " << center.get_atom_type() << std::endl;
        }

        for (auto neigh : center) {
          double distance{
              (center.get_position() - neigh.get_position()).norm()};

          indices_.push_back(neigh.get_atom_tag());
          distances_.push_back(distance);
          auto dir_vec{adaptor_strict->get_direction_vector(neigh)};
          std::array<double, 3> bb{{dir_vec(0), dir_vec(1), dir_vec(2)}};
          dir_vecs_.push_back(bb);

          if (verbose) {
            std::cout << "strict neigh out " << neigh.get_index();
            std::cout << " " << neigh.get_atom_tag() << "\t ";

            for (int ii{0}; ii < 3; ++ii) {
              std::cout << neigh.get_position()[ii] << ", ";
            }
            std::cout << "\t dist=" << distance;
            std::cout << "\t " << neigh.get_atom_type() << std::endl;
          }
        }

        if (verbose) {
          std::cout << "Number of Neighbours: " << indices_.size() << std::endl;
        }

        neigh_ids_strict.push_back(indices_);
        neigh_dist_strict.push_back(distances_);
        neigh_dir_vec_strict.push_back(dir_vecs_);
        // if (icenter > 1) break;
      }

      BOOST_CHECK_EQUAL(neigh_ids.size(), neigh_ids_strict.size());

      for (size_t ii{0}; ii < neigh_ids.size(); ++ii) {
        BOOST_CHECK_EQUAL(neigh_ids[ii].size(), neigh_ids_strict[ii].size());

        for (size_t jj{0}; jj < neigh_ids[ii].size(); ++jj) {
          int a0{neigh_ids[ii][jj]};
          int a1{neigh_ids_strict[ii][jj]};
          double d0{neigh_dist[ii][jj]};
          double d1{neigh_dist_strict[ii][jj]};
          BOOST_CHECK_EQUAL(a0, a1);
          BOOST_CHECK_EQUAL(d0, d1);
          for (size_t kk{0}; kk < neigh_dir_vec[ii][jj].size(); ++kk) {
            double dv0{neigh_dir_vec[ii][jj][kk]};
            double dv1{neigh_dir_vec_strict[ii][jj][kk]};
            BOOST_CHECK_EQUAL(dv0, dv1);
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  // using Fixtures_no_center = boost::mpl::list<
  //     MultipleStructureFixture<MultipleStructureManagerNLCCFixtureCenterMask>,
  //     MultipleStructureFixture<
  //         MultipleStructureManagerNLCCStrictFixtureCenterMask>>;
  using Fixtures_no_center = boost::mpl::list<
      MultipleStructureFixture<MultipleStructureManagerNLCCFixtureCenterMask>>;
  /**
   * Test that + get_size and get_nb_clusters are consistent with their tasks
   *           when masking some atoms
   *
   *           + distances of the strict manager are unchanged, their order
   *           might so sorted distances are compared
   *
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_no_center_test, Fix,
                                   Fixtures_no_center, Fix) {
    bool verbose{false};
    auto & managers = Fix::managers;
    int n_it{static_cast<int>(managers.size())};
    for (int i_it{0}; i_it < n_it; i_it += 2) {
      auto & manager = managers[i_it];
      auto & manager_no_center = managers[i_it + 1];
      auto center_atoms_mask = extract_underlying_manager<0>(manager_no_center)
                                   ->get_center_atoms_mask();

      if (not manager->get_consider_ghost_neighbours()) {
        auto natoms = manager->get_size();
        auto natoms2 = manager->get_nb_clusters(1);
        BOOST_CHECK_EQUAL(natoms, natoms2);

        auto n_center_atom = center_atoms_mask.count();
        natoms = manager_no_center->get_size();
        natoms2 = manager_no_center->get_nb_clusters(1);
        BOOST_CHECK_EQUAL(n_center_atom, natoms);
      } else {
        auto natoms = manager->get_size_with_ghosts();
        auto natoms2 = manager->get_nb_clusters(1);
        BOOST_CHECK_EQUAL(natoms, natoms2);
      }

      if (verbose) {
        std::cout << "center_atoms_mask: " << center_atoms_mask.transpose()
                  << std::endl;
      }
      std::vector<std::vector<double>> distances_ref{};
      std::vector<std::vector<double>> distances{};

      size_t i_center{0};
      if (verbose) {
        std::cout << "Center index: ";
      }
      for (auto center : manager) {
        if (center_atoms_mask(i_center)) {
          distances_ref.emplace_back();
          if (verbose) {
            std::cout << extract_underlying_manager<0>(manager)->get_atom_index(
                             center)
                      << ", ";
          }
          for (auto neigh : center) {
            auto dist{(neigh.get_position() - center.get_position()).norm()};
            distances_ref.back().push_back(dist);
          }
        }
        i_center++;
      }
      if (verbose) {
        std::cout << std::endl;
        std::cout << "Center index center_mask: ";
      }
      for (auto center : manager_no_center) {
        distances.emplace_back();
        if (verbose) {
          std::cout << extract_underlying_manager<0>(manager_no_center)
                           ->get_atom_index(center)
                    << ", ";
        }
        for (auto neigh : center) {
          auto dist{(neigh.get_position() - center.get_position()).norm()};
          distances.back().push_back(dist);
        }
      }
      if (verbose) {
        std::cout << std::endl;
      }
      i_center = 0;
      for (; i_center < manager_no_center->size(); ++i_center) {
        std::sort(distances_ref[i_center].begin(),
                  distances_ref[i_center].end());
        std::sort(distances[i_center].begin(), distances[i_center].end());

        BOOST_CHECK_EQUAL_COLLECTIONS(
            distances_ref[i_center].begin(), distances_ref[i_center].end(),
            distances[i_center].begin(), distances[i_center].end());

        if (verbose) {
          std::cout << "Center: " << i_center << std::endl;
          std::cout << "sizes: " << distances_ref[i_center].size() << ", "
                    << distances[i_center].size() << std::endl;
          for (size_t i_d{0}; i_d < distances[i_center].size(); i_d++) {
            std::cout << std::abs(distances_ref[i_center][i_d] -
                                  distances[i_center][i_d])
                      << "\t" << distances_ref[i_center][i_d] << "\t"
                      << distances[i_center][i_d] << std::endl;
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test that the same hcp structure defined in two ways (conventional and
   * reduced) have the same number of neighbors.
   */
  BOOST_FIXTURE_TEST_CASE(strict_test_hcp, ManagerFixtureTwoHcp) {
    /*
     * Note: since the cell vectors are different, it is possible that one of
     * the two atoms is repeated into a different cell due to periodicity. This
     * leads to a difference in number of neighbours. Therefore the strict
     * cutoff is check to ensure the exact same number of neighbours.
     */
    constexpr bool verbose{false};

    if (verbose) {
      std::cout << "HCP test " << cutoff << std::endl;
    }
    int mult = 3;

    for (auto i{1}; i < mult; ++i) {
      double cutoff_tmp = i * cutoff;

      std::vector<int> neighbours_per_atom1{};
      std::vector<int> neighbours_per_atom2{};

      neighbours_per_atom1.resize(0);
      neighbours_per_atom1.resize(0);

      if (verbose) {
        std::cout << "hcp test cutoff = " << cutoff_tmp << std::endl;
      }

      auto pair_manager1{
          make_adapted_manager<AdaptorNeighbourList>(manager_1, cutoff_tmp)};
      auto adaptor_strict1{
          make_adapted_manager<AdaptorStrict>(pair_manager1, cutoff_tmp)};
      adaptor_strict1->update();

      if (verbose) {
        std::cout << "Setting up strict manager 1 " << std::endl;
      }

      auto pair_manager2{
          make_adapted_manager<AdaptorNeighbourList>(manager_2, cutoff_tmp)};
      auto adaptor_strict2{
          make_adapted_manager<AdaptorStrict>(pair_manager2, cutoff_tmp)};
      adaptor_strict2->update();

      if (verbose) {
        std::cout << "Setting up strict manager 2 " << std::endl;
      }

      for (auto atom : adaptor_strict1) {
        neighbours_per_atom1.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "1 pair " << atom.back() << " " << pair.back()
                      << std::endl;
          }
          adaptor_strict1->get_distance(pair);
          neighbours_per_atom1.back()++;
        }
      }

      for (auto atom : adaptor_strict2) {
        neighbours_per_atom2.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "2 pair " << atom.back() << " " << pair.back()
                      << std::endl;
          }

          neighbours_per_atom2.back()++;
        }
      }

      BOOST_CHECK_EQUAL_COLLECTIONS(
          neighbours_per_atom1.begin(), neighbours_per_atom1.end(),
          neighbours_per_atom2.begin(), neighbours_per_atom2.end());

      for (auto i{0}; i < natoms; ++i) {
        if (verbose) {
          std::cout << "neigh1/neigh2: i " << i << " "
                    << neighbours_per_atom1[i] << "/" << neighbours_per_atom2[i]
                    << std::endl;
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test that the same fcc structure defined in two ways (conventional and
   * reduced) have the same number of neighbors.
   */
  BOOST_FIXTURE_TEST_CASE(neighbourlist_test_fcc, ManagerFixtureTwoFcc) {
    constexpr bool verbose{false};

    if (verbose) {
      std::cout << "FCC test " << std::endl;
    }
    int mult = 3;

    for (auto i{1}; i < mult; ++i) {
      auto cutoff_tmp = i * 0.5 + cutoff;

      std::vector<int> neighbours_per_atom1{};
      std::vector<int> neighbours_per_atom2{};

      neighbours_per_atom1.resize(0);
      neighbours_per_atom1.resize(0);

      if (verbose) {
        std::cout << "fcc cutoff " << cutoff_tmp << std::endl;
      }

      if (verbose) {
        std::cout << "Setting up strict manager 1 " << std::endl;
      }
      auto pair_manager1{
          make_adapted_manager<AdaptorNeighbourList>(manager_1, cutoff_tmp)};
      auto adaptor_strict1{
          make_adapted_manager<AdaptorStrict>(pair_manager1, cutoff_tmp)};
      adaptor_strict1->update();

      if (verbose) {
        std::cout << "Setting up strict manager 2 " << std::endl;
      }
      auto pair_manager2{
          make_adapted_manager<AdaptorNeighbourList>(manager_2, cutoff_tmp)};
      auto adaptor_strict2{
          make_adapted_manager<AdaptorStrict>(pair_manager2, cutoff_tmp)};
      adaptor_strict2->update();

      for (auto atom : adaptor_strict1) {
        neighbours_per_atom1.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "1 pair " << atom.back() << " " << pair.back()
                      << std::endl;
          }
          double dist = {(atom.get_position() - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom1.back()++;
          }
        }
      }

      for (auto atom : adaptor_strict2) {
        neighbours_per_atom2.push_back(0);
        for (auto pair : atom) {
          if (verbose) {
            std::cout << "2 pair " << atom.back() << " " << pair.back()
                      << std::endl;
          }
          double dist = {(atom.get_position() - pair.get_position()).norm()};
          if (dist < cutoff_tmp) {
            neighbours_per_atom2.back()++;
          }
        }
      }

      /*
       * only the first index atom can be checked, since the cell with only one
       * atom does not allow for comparison with other atom's number of
       * neighbours
       */
      BOOST_CHECK_EQUAL(neighbours_per_atom1[0], neighbours_per_atom2[0]);
      if (verbose) {
        std::cout << "neigh1/neigh2: " << neighbours_per_atom1[0] << "/"
                  << neighbours_per_atom2[0] << std::endl;
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
