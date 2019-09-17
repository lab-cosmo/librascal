/**
 * file   test_properties.cc
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  tests for cluster-related properties
 *
 * @section LICENSE
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

#include "test_properties.hh"

namespace rascal {
  // TODO(felix) TODO(alex) test dynamic sized Property completely
  BOOST_AUTO_TEST_SUITE(Property_tests);

  using atom_vector_property_fixtures_with_ghosts =
      OrderOnePropertyBoostList::type_with_ghosts;
  using atom_vector_property_fixtures_without_ghosts =
      OrderOnePropertyBoostList::type_without_ghosts;
  using atom_vector_property_fixtures = OrderOnePropertyBoostList::type;
  using pair_property_fixtures = OrderTwoPropertyBoostList::type;
  using triple_property_fixtures = OrderThreePropertyBoostList::type;

  /* ---------------------------------------------------------------------- */

  BOOST_FIXTURE_TEST_CASE(
      atom_order_one_constructor_tests,
      AtomPropertyFixture<StructureManagerCentersStackFixture>) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << manager->get_name();
      std::cout << " starts now." << std::endl;
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_constructor_tests, Fix,
                                   atom_vector_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(pair_constructor_tests, Fix,
                                   pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(triple_constructor_tests, Fix,
                                   triple_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_FIXTURE_TEST_CASE(
      fill_atom_vector_property_order_one_test,
      AtomPropertyFixture<StructureManagerCentersStackFixture>) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << manager->get_name();
      std::cout << ", manager size " << manager->get_size();
      std::cout << ", manager size with ghosts "
                << manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }

    atom_vector_property.resize(manager->get_consider_ghost_neighbours());
    atom_dynamic_vector_property.resize(
        manager->get_consider_ghost_neighbours());
    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << atom_vector_property.size();
      std::cout << std::endl;
    }
    for (auto atom : manager) {
      if (verbose) {
        std::cout << ">> Atom with tag " << atom.get_atom_tag();
        std::cout << std::endl;
      }
      atom_vector_property[atom] = atom.get_position();
      atom_dynamic_vector_property[atom] = atom.get_position();
    }

    for (auto atom : manager) {
      auto error = (atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      auto error_dynamic =
          (atom_dynamic_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error_dynamic, tol * 100);
    }

    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_sequence_test, Fix,
                                   atom_vector_property_fixtures, Fix) {
    Fix::atom_scalar_property.fill_sequence(
        Fix::manager->get_consider_ghost_neighbours());
    size_t counter{0};
    for (auto atom : Fix::manager->with_ghosts()) {
      BOOST_CHECK_EQUAL(Fix::atom_scalar_property[atom], counter);
      counter++;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_atom_vector_property_test, Fix,
                                   atom_vector_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }

    Fix::atom_scalar_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::atom_vector_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::atom_dynamic_vector_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::atom_dynamic_scalar_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::sparse_atom_scalar_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << std::endl;
    }
    size_t counter{0};
    for (auto atom : Fix::manager->with_ghosts()) {
      if (verbose) {
        std::cout << ">> Atom with tag " << atom.get_atom_tag();
        std::cout << std::endl;
      }
      Fix::atom_vector_property[atom] = atom.get_position();
      Fix::atom_dynamic_vector_property[atom] = atom.get_position();
      Fix::atom_scalar_property[atom] = counter;
      Fix::atom_dynamic_scalar_property[atom] << counter;
      // #TODO(felix) test scalar prop when merged with latest version of
      // the sparse property
      // Fix::sparse_atom_scalar_property[atom] << counter; // DOES NOT WORK
      counter++;
    }

    counter = 0;
    Eigen::MatrixXd eigen_counter(1, 1);
    for (auto atom : Fix::manager->with_ghosts()) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      auto error_dynamic =
          (Fix::atom_dynamic_vector_property[atom] - atom.get_position())
              .norm();
      BOOST_CHECK_LE(error_dynamic, tol * 100);
      BOOST_CHECK_EQUAL(Fix::atom_scalar_property[atom], counter);
      eigen_counter << counter;
      BOOST_CHECK_EQUAL(Fix::atom_dynamic_scalar_property[atom], eigen_counter);
      // BOOST_CHECK_EQUAL(Fix::atom_scalar_property[atom], eigen_counter);
      counter++;
    }

    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  /* ---------------------------------------------------------------------- */
  /**
   * checks if the properties associated with atoms and pairs can be filled
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_pair_property_test, Fix,
                                   pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }

    Fix::atom_vector_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::pair_property.resize();
    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << std::endl;
    }
    int pair_property_counter{};
    for (auto atom : Fix::manager->with_ghosts()) {
      Fix::atom_vector_property[atom] = atom.get_position();
      for (auto pair : atom) {
        Fix::pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom : Fix::manager->with_ghosts()) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(Fix::pair_property[pair], ++pair_property_counter);
      }
    }

    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_triple_property_test, Fix,
                                   triple_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " and consider_ghost_neighbours=";
      std::cout << Fix::manager->get_consider_ghost_neighbours();
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }

    Fix::atom_vector_property.resize(
        Fix::manager->get_consider_ghost_neighbours());
    Fix::pair_property.resize();
    Fix::triple_property.resize();
    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << ", pair_property size ";
      std::cout << Fix::pair_property.size();
      std::cout << ", triple_property size ";
      std::cout << Fix::triple_property.size();
      std::cout << std::endl;
    }
    int pair_property_counter{};
    int triple_property_counter{};
    for (auto atom : Fix::manager->with_ghosts()) {
      Fix::atom_vector_property[atom] = atom.get_position();
      for (auto pair : atom) {
        Fix::pair_property[pair] = ++pair_property_counter;
        for (auto triple : pair) {
          Fix::triple_property[triple] = ++triple_property_counter;
        }
      }
    }

    pair_property_counter = 0;
    triple_property_counter = 0;
    for (auto atom : Fix::manager->with_ghosts()) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(Fix::pair_property[pair], ++pair_property_counter);
        for (auto triple : pair) {
          BOOST_CHECK_EQUAL(Fix::triple_property[triple],
                            ++triple_property_counter);
        }
      }
    }

    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * If consider_ghost_neighbours is true the atoms index should
   * correspond to the cluster index of order 1 when StructureManagerCenters is
   * used as  root implementation and no filtering on order 1 has been done.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_vector_property_fixtures_tests, Fix,
                                   atom_vector_property_fixtures_with_ghosts,
                                   Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    size_t cluster_index{0};
    for (auto atom : Fix::manager) {
      if (verbose) {
        std::cout << ">> Atom index " << atom.get_atom_tag();
        std::cout << ", ClusterIndex should be " << cluster_index;
        std::cout << " and is ";
        std::cout << Fix::manager->get_atom_index(atom.get_atom_tag());
        std::cout << "." << std::endl;
      }
      BOOST_CHECK_EQUAL(Fix::manager->get_atom_index(atom.get_atom_tag()),
                        cluster_index);
      cluster_index++;
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Access of atom property with pair.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_vector_property_access_with_pair_tests,
                                   Fix, pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " and consider_ghost_neighbours=";
      std::cout << Fix::manager->get_consider_ghost_neighbours();
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }
    // initalize the positions
    Fix::atom_scalar_property.resize(false);
    if (verbose) {
      std::cout
          << ">> Property for consider_ghost_atoms=false resized to size ";
      std::cout << Fix::atom_scalar_property.size();
      std::cout << std::endl;
    }
    // Fix::atom_dynamic_scalar_property.resize();
    for (auto atom : Fix::manager) {
      if (verbose) {
        std::cout << ">> Property for atom with tag ";
        std::cout << atom.get_atom_tag();
        std::cout << " is initialized.";
        std::cout << std::endl;
      }
      Fix::atom_scalar_property[atom] = 0;
      // Fix::atom_dynamic_scalar_property[atom] = 0;
      // Fix::sparse_atom_scalar_property[atom] = 0;
    }
    std::vector<size_t> counters{};
    size_t nb_central_atoms = Fix::manager->get_size();
    counters.reserve(nb_central_atoms);
    for (size_t i{0}; i < counters.capacity(); i++) {
      counters.push_back(0);
    }
    if (verbose) {
      std::cout << ">> Counters initialized with size ";
      std::cout << counters.size() << std::endl;
    }
    // add the position to the atom and count how often this happens
    for (auto atom : Fix::manager->with_ghosts()) {
      for (auto pair : atom) {
        if (verbose) {
          std::cout << ">> Atom with tag ";
          std::cout << pair.get_internal_neighbour_atom_tag();
          std::cout << " corresponds to central atom in cell with atom index ";
          std::cout << Fix::manager->get_atom_index(atom.get_atom_tag());
          std::cout << std::endl;
        }
        Fix::atom_scalar_property[pair]++;
        // Fix::atom_dynamic_scalar_property[pair]++;
        // Fix::sparse_atom_scalar_property[pair]++;
        counters.at(Fix::manager->get_atom_index(
            pair.get_internal_neighbour_atom_tag()))++;
      }
    }
    for (auto atom : Fix::manager) {
      size_t counters_at_cluster_index =
          counters.at(Fix::manager->get_atom_index(atom.get_atom_tag()));
      BOOST_CHECK_EQUAL(Fix::atom_scalar_property[atom],
                        counters_at_cluster_index);
      // BOOST_CHECK_EQUAL(Fix::atom_dynamic_scalar_property[atom],
      //    counters_at_cluster_index);
      // BOOST_CHECK_EQUAL(Fix::sparse_atom_scalar_property[atom],
      //    counters_at_cluster_index);
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(order_three_constructor_tests, Fix,
                                   triple_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
      std::cout << " finished." << std::endl;
    }
  }
  /**
   * Access of atom property with triple.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(
      atom_vector_property_access_with_triple_tests, Fix,
      triple_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " and consider_ghost_neighbours=";
      std::cout << Fix::manager->get_consider_ghost_neighbours();
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }
    Fix::atom_scalar_property.resize(false);
    // initalize the positions
    if (verbose) {
      std::cout << ">> atom_vector_property resized to size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << std::endl;
    }
    for (auto atom : Fix::manager) {
      if (verbose) {
        std::cout << ">> Atom tag " << atom.get_atom_tag();
        std::cout << std::endl;
      }
      Fix::atom_scalar_property[atom] = 0;
    }
    std::vector<size_t> counters{};
    size_t nb_central_atoms = Fix::manager->get_size();
    counters.reserve(nb_central_atoms);
    for (size_t i{0}; i < counters.capacity(); i++) {
      counters.push_back(0);
    }

    // add the position to the atom and count how often this happens
    for (auto atom : Fix::manager->with_ghosts()) {
      for (auto pair : atom) {
        for (auto triple : pair) {
          if (verbose) {
            std::cout << ">> Atom with tag "
                      << triple.get_internal_neighbour_atom_tag();
            std::cout << " and cluster index "
                      << Fix::manager->get_atom_index(
                             triple.get_internal_neighbour_atom_tag());
            std::cout << std::endl;
          }
          Fix::atom_scalar_property[triple]++;
          counters.at(Fix::manager->get_atom_index(
              triple.get_internal_neighbour_atom_tag()))++;
        }
      }
    }
    for (auto atom : Fix::manager) {
      if (verbose) {
        std::cout << ">> atom.get_atom_tag() is " << atom.get_atom_tag()
                  << std::endl;
        std::cout << ">> manager->get_atom_index(atom.get_atom_tag()) is "
                  << Fix::manager->get_atom_index(atom.get_atom_tag())
                  << std::endl;
        std::cout << ">> atom_scalar_property[atom] is "
                  << Fix::atom_scalar_property[atom] << std::endl;
        std::cout
            << ">> counters.at(manager->get_atom_index(atom.get_atom_tag())) "
               "is "
            << counters.at(Fix::manager->get_atom_index(atom.get_atom_tag()))
            << std::endl;
      }
      BOOST_CHECK_EQUAL(
          Fix::atom_scalar_property[atom],
          counters.at(Fix::manager->get_atom_index(atom.get_atom_tag())));
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  /* ---------------------------------------------------------------------- */
  /**
   * The access of an order one property with the atom itself
   * and the pair with the atom as neighbour should be the same.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_test_simple_order_one_property, Fix,
                                   pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    Fix::pair_property.resize();
    Fix::atom_vector_property.resize();
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
    }

    for (auto atom : Fix::manager) {
      for (auto pair : atom) {
        if (atom.back() == pair.back()) {
          auto error = (Fix::atom_vector_property[atom] -
                        Fix::atom_vector_property[pair])
                           .norm();
          BOOST_CHECK_LE(error, tol * 100);
        }
      }
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test, if metadata can be assigned to properties
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(meta_data_test, Fix,
                                   atom_vector_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    auto atom_metadata = Fix::atom_vector_property.get_metadata();
    auto dynamic_metadata =
        Fix::atom_dynamic_vector_unit_property.get_metadata();
    auto dynamic_metadata2 = Fix::atom_dynamic_vector_property.get_metadata();

    BOOST_CHECK_EQUAL(atom_metadata, Fix::atom_vector_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata,
                      Fix::atom_dynamic_vector_unit_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata2,
                      Fix::atom_dynamic_vector_property_metadata);
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  /*
   * test filling statically and dynamically sized properties with actual data
   * and comparing if retrieval of those is consistent with the data that was
   * put in
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_test_complex, Fix,
                                   pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    Fix::pair_property.resize();
    Fix::atom_vector_property.resize();
    Fix::atom_dynamic_vector_unit_property.resize();
    Fix::atom_dynamic_vector_property.resize();

    BOOST_CHECK_THROW(Fix::AtomVectorProperty_t ::check_compatibility(
                          Fix::atom_dynamic_vector_unit_property),
                      std::runtime_error);

    BOOST_CHECK_NO_THROW(Fix::AtomVectorProperty_t ::check_compatibility(
        Fix::atom_vector_property));

    int pair_property_counter{};
    size_t counter{};
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
      Fix::atom_dynamic_vector_property[atom] = atom.get_position();

      Fix::atom_dynamic_vector_unit_property[atom] << counter++, counter,
          counter;
      for (auto pair : atom) {
        Fix::pair_property[pair] = ++pair_property_counter;
      }
    }
    auto & base_property{
        static_cast<PropertyBase &>(Fix::atom_dynamic_vector_property)};
    auto & FakeSizedProperty{
        static_cast<typename Fix::AtomVectorProperty_t &>(base_property)};

    pair_property_counter = 0;
    counter = 0;
    for (auto atom : Fix::manager) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      Eigen::Matrix<size_t, Fix::DynSize(), Eigen::Dynamic> tmp(Fix::DynSize(),
                                                                1);
      tmp << counter++, counter, counter;

      auto ierror{(tmp - Fix::atom_dynamic_vector_unit_property[atom]).norm()};
      BOOST_CHECK_EQUAL(ierror, 0);

      error = (Fix::atom_vector_property[atom] -
               Fix::atom_dynamic_vector_property[atom])
                  .norm();
      BOOST_CHECK_LE(error, tol * 100);
      error =
          (Fix::atom_vector_property[atom] - FakeSizedProperty[atom]).norm();
      BOOST_CHECK_LE(error, tol * 100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(Fix::pair_property[pair], ++pair_property_counter);
      }
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  /*
   * test for retrieval of information from property: is it the same that was
   * put in?
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(compute_distances, Fix,
                                   pair_property_fixtures, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    Fix::pair_property.resize();

    for (auto atom : Fix::manager) {
      for (auto pair : atom) {
        Fix::pair_property[pair] =
            (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom : Fix::manager) {
      for (auto pair : atom) {
        auto dist{(atom.get_position() - pair.get_position()).norm()};
        auto error{Fix::pair_property[pair] - dist};
        BOOST_CHECK_LE(error, tol / 100);
      }
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

  /* ---------------------------------------------------------------------- */

  /*
   * A fixture for testing partially sparse properties.
   * TODO(felix) use MultipleStructureManagerCentersFixture instead of NL
   */
  template <size_t Order>
  struct BlockSparsePropertyFixture
      : public MultipleStructureFixture<MultipleStructureManagerNLFixture> {
    using Parent = MultipleStructureFixture<MultipleStructureManagerNLFixture>;
    using ManagerTypeList_t = typename Parent::ManagerTypeHolder_t::type_list;
    using Manager_t = typename Parent::ManagerTypeHolder_t::type;
    using Key_t = std::vector<int>;
    using BlockSparseProperty_t =
        BlockSparseProperty<double, Order, 0, Manager_t, Key_t>;
    using Matrix_t = typename BlockSparseProperty_t::Matrix_t;
    using InputData_t = typename BlockSparseProperty_t::InputData_t;
    using TestData_t = std::vector<InputData_t>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string sparse_features_desc{"some atom centered sparse features"};

    BlockSparsePropertyFixture() : Parent{} {
      std::random_device rd;
      std::mt19937 gen{rd()};
      auto size_dist{std::uniform_int_distribution<size_t>(1, 10)};
      auto key_dist{std::uniform_int_distribution<int>(1, 100)};
      // size_t i_manager{0};
      for (auto & manager : managers) {
        sparse_features.emplace_back(*manager, sparse_features_desc);
        this->keys_list.emplace_back();
        TestData_t test_data{};
        for (size_t i{0}; i < manager->get_size(); ++i) {
          // set up random unique keys
          auto set_max_size{size_dist(gen)};
          std::set<Key_t> keys{};
          for (size_t ii{0}; ii < set_max_size; ii++) {
            keys.emplace(key_dist(gen));
          }

          // set up the data to fill the property later
          InputData_t datas{};
          // resize and set to 0
          datas.resize(keys, n_row, n_col, 0);
          for (auto & key : keys) {
            auto data = Matrix_t::Random(n_row, n_col);
            datas[key] += data;
          }
          this->keys_list.back().push_back(keys);
          test_data.push_back(datas);
        }
        this->test_datas.push_back(test_data);
      }
    }

    int n_row{21};
    int n_col{8};

    std::vector<std::vector<std::set<Key_t>>> keys_list{};
    std::vector<TestData_t> test_datas{};
    std::vector<BlockSparseProperty_t> sparse_features{};
  };

  BOOST_AUTO_TEST_SUITE(Property_partially_sparse_tests);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test, BlockSparsePropertyFixture<1>) {}

  /* ---------------------------------------------------------------------- */
  /*
   * checks if the partially sparse properties associated with centers can be
   * filled and that the data can be accessed consistently.
   */

  using Fixtures = boost::mpl::list<BlockSparsePropertyFixture<1>,
                                    // the order == 2 case test the access of
                                    // the property of order 2 can be properly
                                    // accessed by a clusterRef of order 1
                                    BlockSparsePropertyFixture<2>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_test_simple, Fix, Fixtures, Fix) {
    auto & managers = Fix::managers;
    auto & keys_list = Fix::keys_list;
    auto & sparse_features = Fix::sparse_features;
    auto & n_row = Fix::n_row;
    auto & n_col = Fix::n_col;
    auto & test_datas = Fix::test_datas;

    bool verbose{false};
    // fill the property structures
    auto i_manager{0};
    for (auto & manager : managers) {
      auto & keys{keys_list[i_manager]};
      auto i_center{0};
      sparse_features[i_manager].set_shape(this->n_row, this->n_col);
      sparse_features[i_manager].resize();
      for (auto center : manager) {
        auto && sparse_features_center{sparse_features[i_manager][center]};
        sparse_features_center.resize(keys[i_center], n_row, n_col, 0);
        for (auto & key : keys[i_center]) {
          sparse_features_center[key] = test_datas[i_manager][i_center][key];
        }
        i_center++;
      }
      i_manager++;
    }

    i_manager = 0;
    for (auto & manager : managers) {
      if (verbose)
        std::cout << "manager: " << i_manager << std::endl;
      auto i_center{0};
      for (auto center : manager) {
        if (verbose)
          std::cout << "center: " << i_center << std::endl;

        auto data = sparse_features[i_manager].get_dense_row(center);
        size_t key_id{0};
        double error1{0};
        for (auto & key : sparse_features[i_manager].get_keys(center)) {
          auto && value{test_datas[i_manager][i_center][key]};
          Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> test_data(
              value.data(), value.size());
          error1 += (data.col(key_id) - test_data).squaredNorm();
          key_id += 1;
        }
        if (verbose)
          std::cout << "error1: " << error1 << std::endl;
        BOOST_CHECK_LE(std::sqrt(error1), tol * 100);
        for (auto & key : keys_list[i_manager][i_center]) {
          auto error2 = (sparse_features[i_manager](center, key) -
                         test_datas[i_manager][i_center][key])
                            .norm();
          if (verbose)
            std::cout << "error2: " << error2 << std::endl;
          BOOST_CHECK_LE(error2, tol * 100);
        }
        i_center++;
      }
      i_manager++;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test, if metadata can be assigned to properties
   */
  BOOST_FIXTURE_TEST_CASE(meta_data_test, BlockSparsePropertyFixture<1>) {
    for (auto & sparse_feature : sparse_features) {
      auto sparse_feature_metadata = sparse_feature.get_metadata();
      BOOST_CHECK_EQUAL(sparse_feature_metadata, sparse_features_desc);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
