/**
 * file   test_fields.cc
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  tests for cluster-related fields
 *
 * @section LICENSE
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
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
#include "test_neighbourhood.hh"
#include "neighbourhood_managers/field.hh"


namespace rascal {

  BOOST_AUTO_TEST_SUITE (Field_tests);
  template<class ManagerImplementation>
  struct FieldFixture: public ManagerFixture<ManagerImplementation> {
    using Manager_t = ManagerImplementation;
    using PairScalarField_t = Field<Manager_t, double, 2>;
    using PairScalarField_s_t = Field<Manager_t, double, 3>;
    using AtomVectorField_t = Field<Manager_t, double, 1, 3>;


    FieldFixture()
      :ManagerFixture<ManagerImplementation>{}, pair_field{this->manager},pair_field_s{this->manager},
       atom_field{this->manager}
    {}
    
    PairScalarField_t pair_field;
    PairScalarField_s_t pair_field_s;
    AtomVectorField_t atom_field;
  };

  struct FieldFixture_lammps: public ManagerFixture_lammps {
    using Manager_t = typename ManagerFixture_lammps::Manager_t;
    using PairScalarField_t = Field<Manager_t, double, 2>;
    using AtomVectorField_t = Field<Manager_t, double, 1, 3>;


    FieldFixture_lammps()
      :ManagerFixture_lammps{}, pair_field{this->manager},
       atom_field{this->manager}
    {}

    PairScalarField_t pair_field;
    AtomVectorField_t atom_field;
  };

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_cell, FieldFixture<NeighbourhoodManagerCell>) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(fill_test_cell, FieldFixture<NeighbourhoodManagerCell>) {
    pair_field.resize();
    atom_field.resize();
    int pair_field_counter{};
    for (auto atom: manager) {
      atom_field[atom] = atom.get_position();
      for (auto pair: atom) {
        pair_field[pair] = ++pair_field_counter;
      }
    }

    pair_field_counter = 0;
    for (auto atom: manager) {
      auto error = (atom_field[atom] - atom.get_position()).norm();
      BOOST_CHECK_EQUAL(error, 0);
      for (auto pair: atom) {
        BOOST_CHECK_EQUAL(pair_field[pair], ++pair_field_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(compute_distances_cell, FieldFixture<NeighbourhoodManagerCell>) {
    pair_field.resize();

    for (auto atom: manager) {
      for (auto pair: atom) {
        pair_field[pair] = (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom: manager) {
      for (auto pair: atom) {
        auto error = pair_field[pair] - 1;
        BOOST_CHECK_LE(error, 1e-12);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_strict, FieldFixture<NeighbourhoodManagerStrict>) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(fill_test_strict, FieldFixture<NeighbourhoodManagerStrict>) {
    pair_field_s.resize();
    atom_field.resize();
    int pair_field_counter{};
    for (auto atom: manager) {
      atom_field[atom] = atom.get_position();
      for (auto type: atom) {
        for (auto pair: type) {
          pair_field_s[pair] = ++pair_field_counter;
        }
      }
    }

    pair_field_counter = 0;
    for (auto atom: manager) {
      auto error = (atom_field[atom] - atom.get_position()).norm();
      BOOST_CHECK_EQUAL(error, 0);
      for (auto type: atom) {
        for (auto pair: type) {
          BOOST_CHECK_EQUAL(pair_field_s[pair], ++pair_field_counter);
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(compute_distances_strict, FieldFixture<NeighbourhoodManagerStrict>) {
    pair_field.resize();

    for (auto atom: manager) {
      for (auto type: atom) {
        for (auto pair: type) {
          pair_field_s[pair] = (atom.get_position() - pair.get_position()).norm();
        }
      }
    }

    for (auto atom: manager) {
      for (auto type: atom) {
        for (auto pair: type) {
          auto error = pair_field_s[pair] - 1;
          BOOST_CHECK_LE(error, 1e-12);
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test_lammps, FieldFixture_lammps) {}

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(fill_test_lammps, FieldFixture_lammps) {
    pair_field.resize();
    atom_field.resize();
    int pair_field_counter{};
    for (auto atom: manager) {
      atom_field[atom] = atom.get_position();
      for (auto pair: atom) {
        pair_field[pair] = ++pair_field_counter;
      }
    }

    pair_field_counter = 0;
    for (auto atom: manager) {
      auto error = (atom_field[atom] - atom.get_position()).norm();
      BOOST_CHECK_EQUAL(error, 0);
      for (auto pair: atom) {
        BOOST_CHECK_EQUAL(pair_field[pair], ++pair_field_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(compute_distances_lammps, FieldFixture_lammps) {
    pair_field.resize();

    for (auto atom: manager) {
      for (auto pair: atom) {
        pair_field[pair] = (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom: manager) {
      for (auto pair: atom) {
        auto error = pair_field[pair] - 1;
        BOOST_CHECK_LE(error, 1e-12);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END ();

}  // rascal
