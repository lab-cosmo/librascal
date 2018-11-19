/**
 * file test_adaptor.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   01 Nov 2018
 *
 * @brief Common headers for tests related to `Adaptors`
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

#ifndef TEST_ADAPTOR_H
#define TEST_ADAPTOR_H

#include "tests.hh"
#include "test_structure.hh"

namespace rascal {
  /**
   * This file generates fixtures for testing adators, it is based on previously
   * defined fixtures for `NeighbourHoodManager`s in test_structure.hh
   */
  template<class ManagerImplementation>
  struct PairFixtureSimple : public ManagerFixtureFile<ManagerImplementation>
  {
    using Manager_t = ManagerImplementation;

    static_assert(ManagerImplementation::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    using PairManager_t = AdaptorNeighbourList<ManagerImplementation>;

    PairFixtureSimple()
      : ManagerFixtureFile<ManagerImplementation> {},
      pair_manager{this->manager, 1.}
    {
      this->pair_manager.update();
    }

    ~PairFixtureSimple() {}

    PairManager_t pair_manager;
  };

  /* ---------------------------------------------------------------------- */



}  // rascal


#endif /* TEST_NEIGHBOURHOOD_H */
