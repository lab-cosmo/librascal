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
   * defined fixtures for `NeighbourHoodManager`s in test_structure.hh, it is
   * used in checking the building of the neighbour list in an easy
   * configuration based on 9 atoms.
   */
  template<class ManagerImplementation>
  struct PairFixtureSimple
  {
    using Manager_t = ManagerImplementation;

    static_assert(ManagerImplementation::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    using PairManager_t = AdaptorNeighbourList<ManagerImplementation>;

    PairFixtureSimple():
      cutoff{1.}, pair_manager{fixture.manager, this->cutoff}
    {
      this->pair_manager.update();
    }

    ~PairFixtureSimple() = default;

    ManagerFixtureFile<ManagerImplementation> fixture{};
    double cutoff;
    PairManager_t pair_manager;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * PairFixture based on StructureManagerCenters
   */
  struct PairFixtureCenters
  {
    using Manager_t = StructureManagerCenters;
    using PairManager_t = AdaptorNeighbourList<StructureManagerCenters>;

    static_assert(StructureManagerCenters::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    PairFixtureCenters()
      : cutoff{3.5}, pair_manager{this->fixture.manager, this->cutoff}
    {
      this->pair_manager.update();
    }

    ~PairFixtureCenters() {}

    ManagerFixture<StructureManagerCenters> fixture{};

    double cutoff;
    PairManager_t pair_manager;
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  struct PairFixtureStrict
  {
    using AdaptorStrict_t = AdaptorStrict<ManagerImplementation>;

    PairFixtureStrict():
      adaptor_strict{this->fixture.pair_manager, this->fixture.cutoff}
    {}

    ~PairFixtureStrict() = default;

    // TODO: different fixtures?, streamline fixtures to always work with
    // ´manager´ as an iterator
    PairFixture<ManagerImplementation> fixture{};
    AdaptorStrict_t adaptor_strict;
  };

}  // rascal


#endif /* TEST_NEIGHBOURHOOD_H */
