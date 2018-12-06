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
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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
#include "structure_managers/property.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(Property_tests);

  /* ---------------------------------------------------------------------- */
  /**
   * A fixture for testing proterties. It is based on the PairFixture, which is
   * basically a fixture which provides a pair neighbour list based on positions
   * which are initialized in the tests.
   */
  template<class ManagerImplementation>
  struct PropertyFixture
    : public PairFixture<ManagerImplementation> {
    // TODO(musil) make the type not hard coded
    using Manager_t = AdaptorNeighbourList<ManagerImplementation>;

    using PairScalarProperty_t =
      typename Manager_t::template Property_t<double, 2>;
    using AtomVectorProperty_t =
      typename Manager_t::template Property_t<double, 1, 3, 1>;
    using AtomDynamicProperty_t =
      typename Manager_t::template TypedProperty_t<size_t, 1>;
    using AtomDynamicProperty2_t =
      typename Manager_t::template TypedProperty_t<double, 1>;

    constexpr static Dim_t DynSize() {return 3;}

    std::string atom_property_metadata{"positions"};
    std::string dynamic_property_metadata{"arbitrary counters"};
    std::string dynamic_property2_metadata{"distances"};

    PropertyFixture()
      :PairFixture<ManagerImplementation>{},
      pair_property{this->pair_manager},
      atom_property{this->pair_manager, atom_property_metadata},
      dynamic_property{this->pair_manager, DynSize(), 1,
          dynamic_property_metadata},
      dynamic_property2{this->pair_manager, DynSize(), 1,
          dynamic_property2_metadata}
    {}

    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          PropertyFixture<StructureManagerCenters>) {}

  /* ---------------------------------------------------------------------- */
  /**
   * checks if the properties associated with atoms and pairs can be filled
   */
  BOOST_FIXTURE_TEST_CASE(fill_test_simple,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();
    atom_property.resize();
    int pair_property_counter{};
    for (auto atom : pair_manager) {
      atom_property[atom] = atom.get_position();
      for (auto pair : atom) {
        pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom : pair_manager) {
      auto error = (atom_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol*100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(pair_property[pair], ++pair_property_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test, if metadata can be assigned to properties
   */
  BOOST_FIXTURE_TEST_CASE(meta_data_test,
                          PropertyFixture<StructureManagerCenters> ) {
    auto atom_metadata = atom_property.get_metadata();
    auto dynamic_metadata = dynamic_property.get_metadata();
    auto dynamic_metadata2 = dynamic_property2.get_metadata();

    BOOST_CHECK_EQUAL(atom_metadata, atom_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata, dynamic_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata2, dynamic_property2_metadata);
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test filling statically and dynamically sized properties with actual data
   * and comparing if retrieval of those is consistent with the data that was
   * put in
   */
  BOOST_FIXTURE_TEST_CASE(fill_test_complex,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();
    atom_property.resize();
    dynamic_property.resize();
    dynamic_property2.resize();

    BOOST_CHECK_THROW(AtomVectorProperty_t
                      ::check_compatibility(dynamic_property),
                      std::runtime_error);

    BOOST_CHECK_NO_THROW(AtomVectorProperty_t
                         ::check_compatibility(atom_property));

    int pair_property_counter{};
    size_t counter{};
    for (auto atom : pair_manager) {
      atom_property[atom] = atom.get_position();
      dynamic_property2[atom] = atom.get_position();

      dynamic_property[atom] << counter++, counter, counter;
      for (auto pair : atom) {
        pair_property[pair] = ++pair_property_counter;
      }
    }

    auto & FakeSizedProperty{
      AtomVectorProperty_t::check_compatibility(dynamic_property2)};

    pair_property_counter = 0;
    counter = 0;
    for (auto atom : pair_manager) {
      auto error = (atom_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol*100);
      Eigen::Matrix<size_t, DynSize(), Eigen::Dynamic> tmp(DynSize(), 1);
      tmp << counter++, counter, counter;

      auto ierror{(tmp - dynamic_property[atom]).norm()};
      BOOST_CHECK_EQUAL(ierror, 0);

      error = (atom_property[atom] - dynamic_property2[atom]).norm();
      BOOST_CHECK_LE(error, tol*100);
      error = (atom_property[atom] - FakeSizedProperty[atom]).norm();
      BOOST_CHECK_LE(error, tol*100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(pair_property[pair], ++pair_property_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test for distances
   */
  BOOST_FIXTURE_TEST_CASE(compute_distances,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();

    for (auto atom : pair_manager) {
      for (auto pair : atom) {
        pair_property[pair] =
          (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom : pair_manager) {
      for (auto pair : atom) {
        auto dist{(atom.get_position() - pair.get_position()).norm()};
        auto error = pair_property[pair] - dist;
        BOOST_CHECK_LE(error, tol*100);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
