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
 * Copyright  2018 Markus Stricker, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

    using PairManager_t = AdaptorNeighbourList<Manager_t>;

    PairFixtureSimple()
        : cutoff{1.}, pair_manager{make_adapted_manager<AdaptorNeighbourList>(
                          fixture.manager, this->cutoff)} {}

    ~PairFixtureSimple() = default;

    ManagerFixtureFile<Manager_t> fixture{};
    double cutoff;
    std::shared_ptr<PairManager_t> pair_manager;
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
        : cutoff{3.5}, pair_manager{make_adapted_manager<AdaptorNeighbourList>(
                           this->fixture.manager, this->cutoff, true)} {}

    ~PairFixtureCenters() {}

    ManagerFixture<StructureManagerCenters> fixture{};

    double cutoff;
    std::shared_ptr<PairManager_t> pair_manager;
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  struct PairFixtureStrict {
    using AdaptorStrict_t =
        AdaptorStrict<AdaptorNeighbourList<ManagerImplementation>>;

    PairFixtureStrict()
        : adaptor_strict{make_adapted_manager<AdaptorStrict>(
              this->fixture.pair_manager, this->fixture.cutoff)} {
      this->adaptor_strict->update();
    }

    ~PairFixtureStrict() = default;

    // TODO(markus): different fixtures?, streamline fixtures to always work
    // with ´manager´ as an iterator?
    PairFixture<ManagerImplementation> fixture{};
    std::shared_ptr<AdaptorStrict_t> adaptor_strict;
  };

  /* ---------------------------------------------------------------------- */

  template <class ManagerImplementation>
  struct PairFixtureStrictWithGhosts {
    using AdaptorStrict_t =
        AdaptorStrict<AdaptorNeighbourList<ManagerImplementation>>;

    PairFixtureStrictWithGhosts()
        : adaptor_strict{make_adapted_manager<AdaptorStrict>(
              this->fixture.pair_manager, this->fixture.cutoff)} {
      this->adaptor_strict->update();
    }

    ~PairFixtureStrictWithGhosts() = default;

    // TODO(markus): different fixtures?, streamline fixtures to always work
    // with ´manager´ as an iterator?
    PairFixture<ManagerImplementation> fixture{true};
    std::shared_ptr<AdaptorStrict_t> adaptor_strict;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Streamline the test on several structures and cutoffs
   */

  struct MultipleStructureManagerCentersFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters>;
    MultipleStructureManagerCentersFixture() {
      for (auto && filename : this->filenames) {
        json parameters{};
        json structure{{"filename", filename}};
        json adaptors{};
        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;
        this->factory_args.emplace_back(parameters);
      }
    }

    ~MultipleStructureManagerCentersFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};

    json factory_args{};
  };

  struct MultipleStructureManagerNLFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList>;
    MultipleStructureManagerNLFixture() {
      for (auto && filename : this->filenames) {
        for (auto && cutoff : this->cutoffs) {
          for (auto && skin : this->skins) {
            json parameters;
            json structure{{"filename", filename}};
            json adaptors;
            json ad1{
                {"name", "AdaptorNeighbourList"},
                {"initialization_arguments",
                 {{"cutoff", cutoff},
                  {"skin", skin},
                  {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
            adaptors.push_back(ad1);

            parameters["structure"] = structure;
            parameters["adaptors"] = adaptors;

            this->factory_args.emplace_back(parameters);
          }
        }
      }
    }

    ~MultipleStructureManagerNLFixture() = default;

    const bool consider_ghost_neighbours{false};
    const std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};
    const std::vector<double> cutoffs{{1., 2., 3.}};
    const std::vector<double> skins{{0., 0.3}};

    json factory_args{};
  };

  struct MultipleStructureManagerNLRattleFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList>;
    MultipleStructureManagerNLRattleFixture() {
      for (auto && skin : this->skins) {
        json parameters;
        json structure{{"filename", filename}};
        json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                 {"initialization_arguments",
                  {{"cutoff", cutoff},
                   {"skin", skin},
                   {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
        adaptors.push_back(ad1);

        parameters["structure"] = structure;
        parameters["adaptors"] = adaptors;

        this->factory_args.emplace_back(parameters);
      }
    }

    ~MultipleStructureManagerNLRattleFixture() = default;

    const bool consider_ghost_neighbours{false};
    const std::string filename{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
    const double cutoff{3.};
    const std::vector<double> skins{{0., 0.1, 0.3}};

    json factory_args{};
  };

  struct MultipleStructureManagerNLStrictFixture {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;

    MultipleStructureManagerNLStrictFixture() {
      for (auto && filename : this->filenames) {
        for (auto && cutoff : this->cutoffs) {
          for (auto && skin : this->skins) {
            json parameters;
            json structure{{"filename", filename}};
            json adaptors;
            json ad1{
                {"name", "AdaptorNeighbourList"},
                {"initialization_arguments",
                 {{"cutoff", cutoff},
                  {"skin", skin},
                  {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
            json ad2{{"name", "AdaptorStrict"},
                     {"initialization_arguments", {{"cutoff", cutoff}}}};
            adaptors.emplace_back(ad1);
            adaptors.emplace_back(ad2);

            parameters["structure"] = structure;
            parameters["adaptors"] = adaptors;

            this->factory_args.emplace_back(parameters);
          }
        }
      }
    }

    ~MultipleStructureManagerNLStrictFixture() = default;
    const bool consider_ghost_neighbours{false};
    const std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};
    const std::vector<double> cutoffs{{2., 3.}};
    const std::vector<double> skins{{0., 0.3}};

    json factory_args{};
  };

  template <class BaseFixture>
  struct MultipleStructureFixture : BaseFixture {
    using Parent = BaseFixture;
    using ManagerTypeList_t = typename Parent::ManagerTypeHolder_t::type_list;
    using Manager_t = typename Parent::ManagerTypeHolder_t::type;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    MultipleStructureFixture() : Parent{} {
      for (auto & factory_arg : this->factory_args) {
        auto manager{make_structure_manager_stack_with_hypers_and_typeholder<
            ManagerTypeList_t>::apply(factory_arg["structure"],
                                      factory_arg["adaptors"])};
        this->managers.push_back(manager);
      }
    }

    ~MultipleStructureFixture() = default;

    std::vector<ManagerPtr_t> managers{};
  };

}  // namespace rascal

#endif  // TESTS_TEST_ADAPTOR_HH_
