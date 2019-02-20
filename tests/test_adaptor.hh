/**
 * file test_adaptor.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   01 Nov 2018
 *
 * @brief Common headers for tests related to `Adaptors`
 *
 * @section LICENSE
 *
 * Copyright © 2018 Markus Stricker, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TESTS_TEST_ADAPTOR_HH_
#define TESTS_TEST_ADAPTOR_HH_

#include "tests.hh"
#include "test_structure.hh"

namespace rascal {

  /**
   * This file generates fixtures for testing adators, it is based on previously
   * defined fixtures for `NeighbourHoodManager`s in test_structure.hh, it is
   * used in checking the building of the neighbour list in an easy
   * configuration based on 9 atoms.
   */
  template <class ManagerImplementation>
  struct PairFixtureSimple {
    using Manager_t = ManagerImplementation;

    static_assert(ManagerImplementation::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    using PairManager_t = AdaptorNeighbourList<ManagerImplementation>;

    PairFixtureSimple()
        : cutoff{1.}, pair_manager{fixture.manager, this->cutoff} {
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
  struct PairFixtureCenters {
    using Manager_t = StructureManagerCenters;
    using PairManager_t = AdaptorNeighbourList<Manager_t>;

    static_assert(StructureManagerCenters::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    PairFixtureCenters()
        : cutoff{3.5}, pair_manager{this->fixture.manager, this->cutoff, true} {
      this->pair_manager.update();
    }

    ~PairFixtureCenters() {}

    ManagerFixture<StructureManagerCenters> fixture{};

    double cutoff;
    PairManager_t pair_manager;
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  struct PairFixtureStrict {
    using AdaptorStrict_t = AdaptorStrict<ManagerImplementation>;

    PairFixtureStrict()
        : adaptor_strict{this->fixture.pair_manager, this->fixture.cutoff} {}

    ~PairFixtureStrict() = default;

    // TODO(markus): different fixtures?, streamline fixtures to always work
    // with ´manager´ as an iterator?
    PairFixture<ManagerImplementation> fixture{};
    AdaptorStrict_t adaptor_strict;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Streamline the test on several structures and cutoffs
   */

  struct MultipleStructureManagerBaseFixture {
    MultipleStructureManagerBaseFixture() = default;
    ~MultipleStructureManagerBaseFixture() = default;

    std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};
    std::vector<double> cutoffs{{1., 2., 3.}};
  };

  template <class StructureManager, class BaseFixture>
  struct MultipleStructureManagerCenterFixture : BaseFixture {
    using Parent = BaseFixture;
    using Manager_t = StructureManager;
    using Manager_T_t = AdaptorNeighbourList<Manager_t>;

    MultipleStructureManagerCenterFixture() : Parent{} {
      for (auto filename : this->filenames) {
        this->managers_center.emplace_back();
        this->managers_center.back().update(filename);
      }
    }

    ~MultipleStructureManagerCenterFixture() = default;

    std::list<Manager_t> managers_center{};
  };

  template <class StructureManager, class BaseFixture>
  struct MultipleStructureManagerNLFixture
      : MultipleStructureManagerCenterFixture<StructureManager, BaseFixture> {
    using Parent =
        MultipleStructureManagerCenterFixture<StructureManager, BaseFixture>;
    using Manager_P_t = typename Parent::Manager_t;
    using Manager_t = AdaptorNeighbourList<Manager_P_t>;
    using Manager_T_t = AdaptorStrict<Manager_t>;

    MultipleStructureManagerNLFixture() : Parent{} {
      for (Manager_P_t & manager : this->managers_center) {
        for (double cutoff : this->cutoffs) {
          this->managers_pair.emplace_back(manager, cutoff);
          this->managers_pair.back().update();
        }
      }
    }

    ~MultipleStructureManagerNLFixture() = default;

    std::list<Manager_t> managers_pair{};
  };

  template <class StructureManager, class BaseFixture>
  struct MultipleStructureManagerStrictFixture
      : MultipleStructureManagerNLFixture<StructureManager, BaseFixture> {
    using Parent =
        MultipleStructureManagerNLFixture<StructureManager, BaseFixture>;
    using Manager_P_t = typename Parent::Manager_t;
    using Manager_t = AdaptorStrict<Manager_P_t>;
    // using Manager_T_t = AdaptorStrict<Manager_t>;

    MultipleStructureManagerStrictFixture() : Parent{} {
      for (Manager_P_t & manager : this->managers_pair) {
        this->managers_strict.emplace_back(manager, manager.get_cutoff());
        this->managers_strict.back().update();
      }
    }

    ~MultipleStructureManagerStrictFixture() = default;

    std::list<Manager_t> managers_strict{};
  };

}  // namespace rascal

#endif  // TESTS_TEST_ADAPTOR_HH_
