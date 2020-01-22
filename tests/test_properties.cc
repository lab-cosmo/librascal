/**
 * @file   test_properties.cc
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

#include "rascal/structure_managers/property_lookup_keys.hh"

#include <boost/test/unit_test.hpp>

constexpr double TOLERANCE = 1e-10;

namespace rascal {
  // TODO(felix) TODO(alex) test dynamic sized Property completely
  BOOST_AUTO_TEST_SUITE(property_tests);

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

  /* ---------------------------------------------------------------------- */
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

  /* ---------------------------------------------------------------------- */
  /**
   * Checks for the throw of the runtime error in member .back() of
   * `PropertyTyped`.
   */
  BOOST_FIXTURE_TEST_CASE(
      property_typed_runtime_error,
      AtomPropertyFixture<StructureManagerCentersStackFixture>) {
    BOOST_CHECK_THROW(atom_dynamic_vector_property.back(), std::runtime_error);
  }

  /* ---------------------------------------------------------------------- */
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

    atom_vector_property.resize();
    atom_dynamic_vector_property.resize();

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
      BOOST_CHECK_LE(error, TOLERANCE);
      auto error_dynamic =
          (atom_dynamic_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error_dynamic, TOLERANCE);
    }

    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(fill_sequence_test, Fix,
                                   atom_vector_property_fixtures, Fix) {
    Fix::atom_scalar_property.fill_sequence();
    size_t counter{0};
    for (auto atom : Fix::manager) {
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

    Fix::atom_scalar_property.resize();
    Fix::atom_vector_property.resize();
    Fix::atom_dynamic_vector_property.resize();
    Fix::atom_dynamic_scalar_property.resize();
    Fix::sparse_atom_scalar_property.resize();

    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << std::endl;
    }
    size_t counter{0};
    for (auto atom : Fix::manager) {
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
    for (auto atom : Fix::manager) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, TOLERANCE);
      auto error_dynamic =
          (Fix::atom_dynamic_vector_property[atom] - atom.get_position())
              .norm();
      BOOST_CHECK_LE(error_dynamic, TOLERANCE);
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

    Fix::atom_vector_property.resize();
    Fix::pair_property.resize();
    if (verbose) {
      std::cout << ">> atom_vector_property size ";
      std::cout << Fix::atom_vector_property.size();
      std::cout << std::endl;
    }
    int pair_property_counter{};
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
      for (auto pair : atom.pairs()) {
        Fix::pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom : Fix::manager) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, TOLERANCE);
      for (auto pair : atom.pairs()) {
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
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }

    Fix::atom_vector_property.resize();
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
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
      for (auto pair : atom.pairs()) {
        Fix::pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom : Fix::manager) {
      auto error =
          (Fix::atom_vector_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, TOLERANCE);
      for (auto pair : atom.pairs()) {
        BOOST_CHECK_EQUAL(Fix::pair_property[pair], ++pair_property_counter);
      }
    }

    int triple_property_counter{};
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
      for (auto triple : atom.triplets()) {
        Fix::triple_property[triple] = ++triple_property_counter;
      }
    }

    triple_property_counter = 0;
    for (auto atom : Fix::manager) {
      for (auto triple : atom.triplets()) {
        BOOST_CHECK_EQUAL(Fix::triple_property[triple],
                          ++triple_property_counter);
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
   * The atoms index should correspond to the cluster index of order 1 when
   * StructureManagerCenters is used as  root implementation and no filtering on
   * order 1 has been done.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_vector_property_fixtures_tests, Fix,
                                   atom_vector_property_fixtures, Fix) {
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
      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }
    // initalize the positions
    Fix::atom_scalar_property.resize();
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
    for (auto atom : Fix::manager) {
      for (auto pair : atom.pairs()) {
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

      std::cout << ", manager size " << Fix::manager->get_size();
      std::cout << ", manager size with ghosts "
                << Fix::manager->get_size_with_ghosts();
      std::cout << " starts now." << std::endl;
    }
    Fix::atom_scalar_property.resize();
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
    for (auto atom : Fix::manager) {
      for (auto triplet : atom.triplets()) {
        if (verbose) {
          std::cout << ">> Atom with tag "
                    << triplet.get_internal_neighbour_atom_tag();
          std::cout << " and cluster index "
                    << Fix::manager->get_atom_index(
                           triplet.get_internal_neighbour_atom_tag());
          std::cout << std::endl;
        }
        Fix::atom_scalar_property[triplet]++;
        counters.at(Fix::manager->get_atom_index(
            triplet.get_internal_neighbour_atom_tag()))++;
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
      for (auto pair : atom.pairs()) {
        if (atom.back() == pair.back()) {
          auto error = (Fix::atom_vector_property[atom] -
                        Fix::atom_vector_property[pair])
                           .norm();
          BOOST_CHECK_LE(error, TOLERANCE);
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

    BOOST_CHECK_THROW(Fix::AtomVectorProperty_t::check_compatibility(
                          Fix::atom_dynamic_vector_unit_property),
                      std::runtime_error);

    BOOST_CHECK_NO_THROW(Fix::AtomVectorProperty_t::check_compatibility(
        Fix::atom_vector_property));

    int pair_property_counter{};
    size_t counter{};
    for (auto atom : Fix::manager) {
      Fix::atom_vector_property[atom] = atom.get_position();
      Fix::atom_dynamic_vector_property[atom] = atom.get_position();

      Fix::atom_dynamic_vector_unit_property[atom] << counter++, counter,
          counter;
      for (auto pair : atom.pairs()) {
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
      BOOST_CHECK_LE(error, TOLERANCE);
      Eigen::Matrix<size_t, Fix::DynSize(), Eigen::Dynamic> tmp(Fix::DynSize(),
                                                                1);
      tmp << counter++, counter, counter;

      auto ierror{(tmp - Fix::atom_dynamic_vector_unit_property[atom]).norm()};
      BOOST_CHECK_EQUAL(ierror, 0);

      error = (Fix::atom_vector_property[atom] -
               Fix::atom_dynamic_vector_property[atom])
                  .norm();
      BOOST_CHECK_LE(error, TOLERANCE);
      error =
          (Fix::atom_vector_property[atom] - FakeSizedProperty[atom]).norm();
      BOOST_CHECK_LE(error, TOLERANCE);
      for (auto pair : atom.pairs()) {
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
      for (auto pair : atom.pairs()) {
        Fix::pair_property[pair] =
            (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom : Fix::manager) {
      for (auto pair : atom.pairs()) {
        auto dist{(atom.get_position() - pair.get_position()).norm()};
        auto error{Fix::pair_property[pair] - dist};
        BOOST_CHECK_LE(error, TOLERANCE);
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
  /** Tests the utility functions in structure manager file including
   * get_property and create property.
   */
  BOOST_AUTO_TEST_SUITE(property_structure_manager_property_tests);

  BOOST_FIXTURE_TEST_CASE(create_property,
                          StructureManagerCentersStackFixture) {
    BOOST_CHECK_NO_THROW(
        this->manager->template create_property<Property_t>("prop"));
    BOOST_CHECK_THROW(
        this->manager->template create_property<Property_t>("prop"),
        std::runtime_error);
  }

  /**
   * Tests if property forwarding works for a stack with two layers
   *
   * ANL -> has not "prop"
   *  | forwards request to lower stack
   *  v
   * SMC has "prop"
   */
  BOOST_FIXTURE_TEST_CASE(property_forwarding_two_layers,
                          ANLWithGhosts_SMC_StackFixture) {
    using Parent_Property_t = Parent::Property_t;
    this->manager->get_previous_manager()
        ->template create_property<Parent_Property_t>("prop");
    BOOST_CHECK_NO_THROW(
        this->manager->template get_property<Property_t>("prop", false););
  }
  /**
   * Tests if property forwarding works for a stack with three layers
   *
   * AS -> has not "prop"
   *  | forwards request to lower stack
   * ANL -> has not "prop"
   *  | forwards request to lower stack
   *  v
   * SMC has "prop"
   */
  BOOST_FIXTURE_TEST_CASE(
      property_forwarding_three_layers,
      AdaptorStrictStackFixture<ANLWithGhosts_SMC_StackFixture>) {
    using Parent_Parent_Property_t = Parent::Parent::Property_t;
    // using Parent_Parent_Property_t = Parent_Parent_t::Property_t;
    this->manager->get_previous_manager()
        ->get_previous_manager()
        ->template create_property<Parent_Parent_Property_t>("prop");
    BOOST_CHECK_NO_THROW(
        this->manager->template get_property<Property_t>("prop", false););
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
        BlockSparseProperty<double, Order, Manager_t, Key_t>;
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

  using Fixtures = boost::mpl::list<BlockSparsePropertyFixture<1>>;

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
        BOOST_CHECK_LE(std::sqrt(error1), TOLERANCE);
        for (auto & key : keys_list[i_manager][i_center]) {
          auto error2 = (sparse_features[i_manager](center, key) -
                         test_datas[i_manager][i_center][key])
                            .norm();
          if (verbose)
            std::cout << "error2: " << error2 << std::endl;
          BOOST_CHECK_LE(error2, TOLERANCE);
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

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE(tests_without_reference_data);

  BOOST_FIXTURE_TEST_CASE(caller_layer_test,
                          ManagerFixture<StructureManagerLammps>) {
    bool verbose{false};

    double cutoff{1.1};  // one pair dist is sqrt{2}, the others are 1
    auto strict{make_adapted_manager<AdaptorStrict>(this->manager, cutoff)};
    strict->update();

    constexpr static int PairOrder{2};
    using PropertyLow_t =
        typename Manager_t::template Property_t<int, PairOrder>;

    auto & low_prop{
        this->manager->template create_property<PropertyLow_t>("low layer")};
    low_prop.resize();

    using PropertyHigh_t =
        typename AdaptorStrict<Manager_t>::template Property_t<int, PairOrder>;

    auto & high_prop{
        strict->template create_property<PropertyHigh_t>("high layer")};
    high_prop.resize();

    int counter{0};
    for (auto && atom : this->manager) {
      auto x{atom.get_position().transpose()};
      for (auto && pair : atom.pairs()) {
        auto y{pair.get_position().transpose()};
        auto dist = (x - y).norm();
        low_prop[pair] = counter;
        counter++;
        if (verbose) {
          std::cout << "low: " << low_prop[pair] << ", dist = " << dist
                    << ", pos a: " << x << ", pos b: " << y << std::endl;
        }
      }
    }

    counter = 0;
    for (auto && atom : *strict) {
      for (auto && pair : atom.pairs()) {
        high_prop[pair] = counter;
        counter++;
        if (verbose) {
          std::cout << "high: " << high_prop[pair]
                    << ", dist = " << strict->get_distance(pair)
                    << ", pos a: " << atom.get_position().transpose()
                    << ", pos b: " << pair.get_position().transpose()
                    << std::endl;
        }
      }
    }

    // strict is supposed to contain pairs 1 and 3
    std::array<int, 2> low_pairs{{1, 3}};
    std::array<int, 2> high_pairs{{0, 1}};
    counter = 0;
    for (auto && atom : *strict) {
      for (auto && pair : atom.pairs()) {
        BOOST_CHECK_EQUAL(low_prop[pair], low_pairs[counter]);
        BOOST_CHECK_EQUAL(high_prop[pair], high_pairs[counter]);
        counter++;
        if (verbose) {
          std::cout << "low: " << low_prop[pair]
                    << ", pos a: " << atom.get_position().transpose()
                    << ", pos b: " << pair.get_position().transpose()
                    << std::endl;
          std::cout << "high: " << high_prop[pair]
                    << ", pos a: " << atom.get_position().transpose()
                    << ", pos b: " << pair.get_position().transpose()
                    << std::endl;
        }
      }
    }

    // strict is supposed to contain pairs 1 and 3
    auto & low_prop_typed{dynamic_cast<TypedProperty<
        int, PairOrder, StructureManager<StructureManagerLammps>> &>(low_prop)};
    auto & high_prop_typed{dynamic_cast<TypedProperty<
        int, PairOrder,
        StructureManager<std::remove_reference_t<decltype(*strict)>>> &>(
        high_prop)};
    counter = 0;
    for (auto && atom : *strict) {
      for (auto && pair : atom.pairs()) {
        BOOST_CHECK_EQUAL(low_prop_typed[pair](0), low_pairs[counter]);
        BOOST_CHECK_EQUAL(high_prop_typed[pair](0), high_pairs[counter]);
        counter++;
        if (verbose) {
          std::cout << "low: " << low_prop_typed[pair]
                    << ", pos a: " << atom.get_position().transpose()
                    << ", pos b: " << pair.get_position().transpose()
                    << std::endl;
          std::cout << "high: " << high_prop_typed[pair]
                    << ", pos a: " << atom.get_position().transpose()
                    << ", pos b: " << pair.get_position().transpose()
                    << std::endl;
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(storing_cluster_ref_keys_for_lookup_test,
                          ManagerFixtureSimple4) {
    auto double_pair_manager{make_adapted_manager<AdaptorNeighbourList>(
        this->manager, this->cutoff)};
    auto loose_pair_manager{
        make_adapted_manager<AdaptorHalfList>(double_pair_manager)};
    auto pair_manager{
        make_adapted_manager<AdaptorStrict>(loose_pair_manager, this->cutoff)};
    using PairManager_t = typename decltype(pair_manager)::element_type;
    auto triplet_manager{
        make_adapted_manager<AdaptorMaxOrder>(pair_manager)};

    triplet_manager->update();
    auto & pair_to_i_atom{
        triplet_manager->template get_neighbours_to_i_atoms<PairOrder>()};
    auto & trip_to_i_atom{
        triplet_manager->template get_neighbours_to_i_atoms<TripletOrder>()};

    constexpr auto AtomLayer{
        decltype(triplet_manager)::element_type::
            template cluster_layer_from_order<AtomOrder>()};
    using AtomClusterRef_t = ClusterRefKey<AtomOrder, AtomLayer>;
    constexpr auto PairLayer{
        decltype(triplet_manager)::element_type::
            template cluster_layer_from_order<PairOrder>()};
    using StoredRefKey_t = ClusterRefKey<PairOrder, PairLayer>;

    using Key_t = std::array<AtomClusterRef_t, PairOrder>;
    std::map<Key_t, StoredRefKey_t> reverse_map{};


    for (auto && atom : pair_manager) {
      for (auto && pair : atom.pairs()) {
        auto && atom_cluster_indices{pair_to_i_atom[pair]};
        auto && i_atom_id{atom_cluster_indices(0)};
        auto && j_atom_id{atom_cluster_indices(1)};

        Key_t key{*triplet_manager->get_iterator_at(i_atom_id),
                  *triplet_manager->get_iterator_at(j_atom_id)};
        std::pair<Key_t, StoredRefKey_t> key_val(
            key, dynamic_cast<StoredRefKey_t &>(pair));

        reverse_map.insert(key_val);
      }
    }

    using Clusters =
        PropertyLookupKeys<StoredRefKey_t, TripletOrder,
                           decltype(triplet_manager)::element_type, 2>;

    auto & clusters{triplet_manager->template create_property<Clusters>(
        "pair_cluster_indices")};

    clusters.resize();

    auto & r_ij_r_ik_r_jk{triplet_manager->template create_property<Property<
        double, TripletOrder, decltype(triplet_manager)::element_type, 3>>(
        "r_ij_r_ik_r_jk")};
    r_ij_r_ik_r_jk.resize();
    auto & distances = triplet_manager.get_distance();
    for (auto && atom : triplet_manager) {
      for (auto && trip : atom.triplets()) {
        auto && atom_cluster_indices{trip_to_i_atom[trip]};
        auto && i_atom_id{atom_cluster_indices(0)};
        auto && j_atom_id{atom_cluster_indices(1)};
        auto && k_atom_id{atom_cluster_indices(2)};
        std::cout << trip.get_atom_tag_list()[0] << ", "
                  << trip.get_atom_tag_list()[1] << ", "
                  << trip.get_atom_tag_list()[2] << std::endl;

        AtomClusterRef_t i_atom{triplet_manager->operator[](i_atom_id)};
        AtomClusterRef_t j_atom{triplet_manager->operator[](j_atom_id)};
        AtomClusterRef_t k_atom{triplet_manager->operator[](k_atom_id)};

        auto & p_ij{reverse_map[{i_atom, j_atom}]};
        auto & p_ik{reverse_map[{i_atom, k_atom}]};

        auto && pairs{clusters[trip]};
        pairs[0] = p_ij;
        pairs[1] = p_ik;

        auto && dists{r_ij_r_ik_r_jk[trip]};
        dists(0) = distances[p_ij];
        dists(1) = distances[p_ik];
        dists(2) = (triplet_manager->position(trip.get_atom_tag_list()[2]) -
                    triplet_manager->position(trip.get_atom_tag_list()[1]))
                       .norm();
      }
    }

    auto & r_ij_r_ik_r_jk_ref{
        triplet_manager->template create_property<Property<
            double, TripletOrder, decltype(triplet_manager)::element_type, 3>>(
            "r_ij_r_ik_r_jk_ref")};

    using Dists_t = Eigen::Matrix<double, 3, 1>;
    r_ij_r_ik_r_jk_ref.push_back(Dists_t::Zero());
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
