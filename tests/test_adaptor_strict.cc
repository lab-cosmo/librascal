/**
 * file   test_adaptor_strict.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
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

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"

#include <vector>

namespace rascal {

  using Fixtures = boost::mpl::list<
      PairFixtureStrict<PairFixtureSimple<StructureManagerCenters>>,
      PairFixtureStrict<PairFixtureCenters>
      // PairFixtureStrict<ManagerFixture<StructureManagerLammps>>
      >;

  BOOST_AUTO_TEST_SUITE(strict_adaptor_test);

  /* ---------------------------------------------------------------------- */
  // TODO(markus): some template parameter does not work with this Fixture list
  // BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Fixtures, Fix) {
  // }

  /* ---------------------------------------------------------------------- */
  /**
   * test strict neighbourhood constructor, this is done here instead
   * of in the Fixture, because of the default constructor is deleted.
   */
  BOOST_FIXTURE_TEST_CASE(constructor_test, PairFixtureCenters) {
    auto adaptor_strict{
        make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Update test
   */
  BOOST_FIXTURE_TEST_CASE(update_test, PairFixtureCenters) {
    auto adaptor_strict{
        make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
    adaptor_strict->update();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Iteration test for strict adaptor. It also checks if the positions and
   * types of the original atoms are accessed correctly.
   */
  BOOST_FIXTURE_TEST_CASE(iterator_test, PairFixtureCenters) {
    auto adaptor_strict{
        make_adapted_manager<AdaptorStrict>(pair_manager, cutoff)};
    adaptor_strict->update();

    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{false};

    for (auto atom : adaptor_strict) {
      auto index{atom.get_global_index()};
      BOOST_CHECK_EQUAL(index, atom_counter);

      auto type{atom.get_atom_type()};
      BOOST_CHECK_EQUAL(type, this->fixture.atom_types[index]);
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
  /* ---------------------------------------------------------------------- */
  /*
   * Compare the strict neighbour list with the linked cell one
   * selecting only the atoms within a cutoff radius
   *
   * ``manager`` is a MaxOrder=1 StructureManager
   */
  BOOST_FIXTURE_TEST_CASE(strict_test,
                          ManagerFixture<StructureManagerCenters>) {
    bool verbose{false};
    int mult = 3;
    // double rc_max{mult * 0.5 + cutoff};
    // auto pair_manager{make_adapted_manager<AdaptorNeighbourList>(manager,
    // rc_max)}; pair_manager->update();

    for (auto i{0}; i < mult; ++i) {
      auto cutoff_tmp{i * 0.5 + cutoff};
      std::vector<std::vector<int>> neigh_ids{};
      std::vector<std::vector<double>> neigh_dist{};
      std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec{};
      std::vector<std::vector<int>> neigh_ids_strict{};
      std::vector<std::vector<double>> neigh_dist_strict{};
      std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec_strict{};

      if (verbose) {
        std::cout << "Setting up strict manager with rc = " << cutoff_tmp
                  << std::endl;
      }
      auto pair_manager{
          make_adapted_manager<AdaptorNeighbourList>(manager, cutoff_tmp)};
      auto adaptor_strict{
          make_adapted_manager<AdaptorStrict>(pair_manager, cutoff_tmp)};
      adaptor_strict->update();

      if (verbose) {
        std::cout << "Setting up comparison list with rc = " << cutoff_tmp
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
          if (distance <= cutoff_tmp) {
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
          double distance{adaptor_strict->get_distance(neigh)};

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
          std::cout << "Number of Neighbourg: " << indices_.size() << std::endl;
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

  using multiple_fixtures = boost::mpl::list<
      MultipleStructureFixture<MultipleStructureManagerNLFixture>>;
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
  BOOST_FIXTURE_TEST_CASE(strict_test_hcp, ManagerFixtureTwoHcp) {
    /*
     * Note: since the cell vectors are different, it is possible that one of
     * the two atoms is repeated into a different cell due to periodicity. This
     * leads to a difference in number of neighbours. Therefore the strict
     * cutoff is check to ensure the exakt same number of neighbours.
     */
    constexpr bool verbose{false};

    if (verbose) {
      std::cout << "HCP test " << cutoff << std::endl;
    }
    int mult = 3;

    for (auto i{1}; i < mult; ++i) {
      auto cutoff_tmp = i * 0.5 + cutoff;

      std::vector<int> neighbours_per_atom1{};
      std::vector<int> neighbours_per_atom2{};

      neighbours_per_atom1.resize(0);
      neighbours_per_atom1.resize(0);

      if (verbose) {
        std::cout << "hcp test cutoff " << cutoff_tmp << std::endl;
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
