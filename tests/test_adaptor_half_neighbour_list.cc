/**
 * file   test_adaptor_half_neighbour_list.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief  tests the implementation of the half neighbourlist adaptor
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

  BOOST_AUTO_TEST_SUITE(half_neighbourlist_adaptor_test);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerLammps>) {

    //! TODO: should this be included in the constructor?
    // I don't think so because if it fails in the Fixture is is harder
    // to undertand where the error comes from
    AdaptorHalfList<StructureManagerLammps> adaptor{manager};
    adaptor.update();
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(iteration_and_distance_half_list,
                          ManagerFixture<StructureManagerLammps>) {

    constexpr bool verbose{false};

    double distance_sum_full{0.};

    int npairs_full{0};
    for (auto atom : manager) {
      for (auto pair : atom) {
        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_full += dist;
        npairs_full++;
      }
    }

    if (verbose) {
      std::cout << "Setting up half neighbourlist manager" << std::endl;
    }

    AdaptorHalfList<StructureManagerLammps> adaptor{manager};
    adaptor.update();

    double distance_sum_half{0.};

    int npairs_half{0};
    for (auto atom : adaptor) {
      if (verbose) std::cout << "type " << atom.get_atom_type() << std::endl;
      for (auto pair : atom) {

        double dist = {(atom.get_position() - pair.get_position()).norm()};
        distance_sum_half += dist;

        dist = {(pair.get_position() - atom.get_position()).norm()};
        distance_sum_half += dist;
        npairs_half++;
      }
    }

    if (verbose) {
      std::cout << "Full/half "
                << npairs_full
                << "/"
                << npairs_half << std::endl;
    }

    auto val{distance_sum_full - distance_sum_half};
    auto relative_error = val * val / (distance_sum_full * distance_sum_full);

    BOOST_CHECK_EQUAL(npairs_full, 4);
    BOOST_CHECK_EQUAL(npairs_half, 2);
    BOOST_CHECK(relative_error < tol * tol);

  }

  /* ---------------------------------------------------------------------- */
  // TODO: add test for going to full neighbour list again

  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
