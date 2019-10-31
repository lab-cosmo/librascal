/**
 * @file   test_atomic_structure.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   3 June 2019
 *
 * @brief test atomic_structure class
 *
 * Copyright © 2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "atomic_structure.hh"

#include <boost/test/unit_test.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(AtomicStructureTests);

  struct AtomicStructureFixture {
    AtomicStructureFixture() {}

    ~AtomicStructureFixture() = default;

    std::string ref_filename1 =
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json";
    std::string ref_filename2 = "reference_data/small_molecule.json";

    bool verbose{false};
  };
  /* ---------------------------------------------------------------------- */
  /**
   * Test the loading of a structure from a json file and the test for
   * identity between structures.
   */
  BOOST_FIXTURE_TEST_CASE(atomic_structure_test, AtomicStructureFixture) {
    AtomicStructure<3> structure1{};
    AtomicStructure<3> structure2{};
    AtomicStructure<3> structure3{};

    // load structure from a json formated file
    structure1.set_structure(ref_filename1);
    structure2.set_structure(ref_filename2);

    // check if identical with itself
    double skin2{0.};
    BOOST_CHECK(structure1.is_similar(structure1, skin2));
    BOOST_CHECK(structure2.is_similar(structure2, skin2));
    BOOST_CHECK(not structure1.is_similar(structure2, skin2));

    skin2 = 0.1 * 0.1;
    structure3.set_structure(structure1);
    structure3.pbc(0) = false;
    BOOST_CHECK(not structure1.is_similar(structure3, skin2));

    structure3.set_structure(structure1);
    structure3.cell(0, 0) = 20;
    BOOST_CHECK(not structure1.is_similar(structure3, skin2));

    structure3.set_structure(structure1);
    structure3.positions(0, 0) += 0.05;
    BOOST_CHECK(structure1.is_similar(structure3, skin2));

    structure3.set_structure(structure1);
    structure3.positions(0, 0) += 0.15;
    BOOST_CHECK(not structure1.is_similar(structure3, skin2));

    structure3.set_structure(structure1);
    structure3.positions(0, 0) += 0.1;
    BOOST_CHECK(not structure1.is_similar(structure3, skin2));
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Test the wrapping of the atoms in a structure
   */
  BOOST_FIXTURE_TEST_CASE(wrap_positions_test, AtomicStructureFixture) {
    AtomicStructure<3> structure1{};
    AtomicStructure<3> structure2{};
    bool verbose{false};

    // load structure from a json formated file
    structure1.set_structure(
        std::string("./reference_data/dummy_structure.json"));
    structure1.wrap();
    structure2.set_structure(
        std::string("./reference_data/dummy_structure_wrapped.json"));

    // check if identical with itself
    double skin2{1e-15};
    BOOST_CHECK(structure1.is_similar(structure2, skin2));
    if (verbose) {
      std::cout << (structure2.positions - structure1.positions).transpose()
                << std::endl;
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
