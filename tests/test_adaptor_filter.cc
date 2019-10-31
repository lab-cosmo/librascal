/**
 * @file   test_adaptor_filter.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   26 Nov 2018
 *
 * @brief  Test the filter adaptor
 *
 * Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/adaptor_filter.hh"
#include "test_structure.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <random>

namespace rascal {

  template <class ManagerImplementation, size_t MaxOrder>
  struct FilterFixture {
    class Filter_t : public AdaptorFilter<ManagerImplementation, MaxOrder> {
      using AdaptorFilter<ManagerImplementation, MaxOrder>::AdaptorFilter;
      void perform_filtering() final{};
    };

    FilterFixture() : manager{fixture.manager} {}

    ~FilterFixture() = default;

    size_t get_MaxOrder() { return MaxOrder; }

    ManagerFixture<ManagerImplementation> fixture{};
    Filter_t manager;
  };

  using Fixtures = boost::mpl::list<FilterFixture<StructureManagerCenters, 1>,
                                    FilterFixture<StructureManagerLammps, 2>>;

  using FixturesMax1 =
      boost::mpl::list<FilterFixture<StructureManagerCenters, 1>>;
  using FixturesMax2 =
      boost::mpl::list<FilterFixture<StructureManagerLammps, 2>>;
  using FixturesMax3 = boost::mpl::list<>;

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE(test_adaptor_filter);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Fixtures, Fix) {}

  /**
   * This test iterates over a Structure manager's atoms, randomly
   * filters them or not using AdaptorFilter, and stores the indices
   * of the retained atoms in a separate backup vector. In a second
   * loop, we check that the content of the AdaptorFilter is identical
   * to the backup vector. The random filtering is meant to catch all
   * corner cases which might have been missed the authors.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(filter_1_test, Fix, FixturesMax1, Fix) {
    std::random_device rd{};
    std::uniform_int_distribution<int> dist(0, 1);
    std::vector<int> atom_tag_list{};

    for (auto atom : Fix::fixture.manager) {
      const bool include(dist(rd));
      if (include) {
        Fix::manager.add_cluster(atom);
        atom_tag_list.push_back(atom.get_atom_tag());
      }
    }

    size_t counter{0};
    for (auto atom : Fix::manager) {
      BOOST_CHECK_EQUAL(atom.get_atom_tag(), atom_tag_list[counter]);
      counter++;
      const auto & pos_a{atom.get_position()};
      const auto & pos_b{
          this->fixture.manager->get_position(atom.get_atom_tag())};
      const auto error{(pos_a - pos_b).norm()};
      BOOST_CHECK_EQUAL(error, 0.);

      const auto & atom_type_a{atom.get_atom_type()};
      const auto & atom_type_b{
          this->fixture.manager->get_atom_type(atom.back())};
      BOOST_CHECK_EQUAL(atom_type_a, atom_type_b);
    }
  }

  /**
   * This test iterates over a Structure manager's pairs, randomly
   * filters them or not using AdaptorFilter, and stores the atom
   * indices of the retained atom pairs in a separate backup
   * vector. In a second loop, we check that the content of the
   * AdaptorFilter is identical to the backup vector. The random
   * filtering is meant to catch all corner cases which might have
   * been missed the authors.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(filter_2_test, Fix, FixturesMax2, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::fixture.manager->get_name();
      std::cout << ", manager size " << Fix::fixture.manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::fixture.manager->get_size_with_ghosts();
      std::cout << ", manager size nb_clusters[Order=2] "
                << Fix::fixture.manager->get_nb_clusters(2);
      std::cout << " starts now." << std::endl;
    }
    std::random_device rd{};
    std::uniform_int_distribution<int> dist(0, 1);
    std::vector<std::array<int, 2>> atom_tag_list{};

    for (auto atom : Fix::fixture.manager) {
      if (verbose) {
        std::cout << ">> Atom with cluster index ";
        std::cout << atom.get_cluster_indices()[0];
        std::cout << " and with atom tag ";
        std::cout << atom.get_atom_tag();
        std::cout << std::endl;
      }
      for (auto pair : atom) {
        if (verbose) {
          std::cout << ">> Pair with cluster index ";
          std::cout << pair.get_cluster_indices()[0];
          std::cout << " and with atom tags ";
          std::cout << pair.get_atom_tag_list()[0];
          std::cout << " ";
          std::cout << pair.get_atom_tag_list()[1];
          std::cout << std::endl;
        }
        const bool include(dist(rd));
        if (include) {
          Fix::manager.add_cluster(pair);
          atom_tag_list.push_back(pair.get_atom_tag_list());
        }
      }
    }

    if (verbose) {
      std::cout << ">> Pairs added" << std::endl;
    }

    size_t counter{0};
    for (auto atom : Fix::manager) {
      if (verbose) {
        std::cout << ">> Atom with cluster index ";
        std::cout << atom.get_cluster_indices()[0];
        std::cout << " and with atom tag ";
        std::cout << atom.get_atom_tag();
        std::cout << std::endl;
      }
      for (auto pair : atom) {
        if (verbose) {
          std::cout << ">> Pair with cluster index ";
          std::cout << pair.get_cluster_indices()[0];
          std::cout << " and with atom tags ";
          std::cout << pair.get_atom_tag_list()[0];
          std::cout << " ";
          std::cout << pair.get_atom_tag_list()[1];
          std::cout << std::endl;
        }
        auto && a{pair.get_atom_tag_list()};
        auto && b{atom_tag_list[counter]};
        BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), b.begin(), b.end());

        const auto & pos_a{pair.get_position()};
        const auto & pos_b{
            this->fixture.manager->get_position(pair.get_atom_tag())};
        const auto error{(pos_a - pos_b).norm()};
        BOOST_CHECK_EQUAL(error, 0.);

        const auto & atom_type_a{pair.get_atom_type()};
        const auto & atom_type_b{
            this->fixture.manager->get_atom_type(pair.back())};
        BOOST_CHECK_EQUAL(atom_type_a, atom_type_b);
        counter++;
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
