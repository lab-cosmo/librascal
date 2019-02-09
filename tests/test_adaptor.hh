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
  template <class ManagerImplementation>
  struct PairFixtureSimple {
    using Manager_t = ManagerImplementation;

    static_assert(ManagerImplementation::traits::MaxOrder == 1,
                  "Lower layer manager has to be a collection of atoms, i.e."
                  " MaxOrder=1");

    using PairManager_t = AdaptorNeighbourList<Manager_t>;

    PairFixtureSimple()
        : cutoff{1.}, pair_manager{make_adapted_manager<AdaptorNeighbourList>(fixture.manager, this->cutoff)} {
      fixture.manager->add_child(this->pair_manager);
      this->pair_manager->update();
    }

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
        : cutoff{3.5}, pair_manager{make_adapted_manager<AdaptorNeighbourList>(this->fixture.manager, this->cutoff, true)} {
      fixture.manager->add_child(this->pair_manager);
      this->pair_manager->update();
    }

    ~PairFixtureCenters() {}

    ManagerFixture<StructureManagerCenters> fixture{};

    double cutoff;
    std::shared_ptr<PairManager_t> pair_manager;
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  struct PairFixtureStrict {
    using AdaptorStrict_t = AdaptorStrict<ManagerImplementation>;

    PairFixtureStrict()
        : adaptor_strict{make_adapted_manager<AdaptorStrict>(this->fixture.pair_manager, this->fixture.cutoff)} {
          fixture.pair_manager->add_child(this->adaptor_strict);
        }

    ~PairFixtureStrict() = default;

    // TODO(markus): different fixtures?, streamline fixtures to always work
    // with ´manager´ as an iterator?
    PairFixture<ManagerImplementation> fixture{};
    std::shared_ptr<AdaptorStrict_t> adaptor_strict;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Streamline the test on several structures and cutoffs
   */

  struct MultipleStructureManagerNLFixture {
    using Factory_t = std::tuple<std::string,std::tuple<double>>;
    using AdaptorTypeHolder_t = AdaptorTypeHolder<StructureManagerCenters, AdaptorNeighbourList>;
    MultipleStructureManagerNLFixture() {
      for (auto&& filename : this->filenames) {
        for (auto&& cutoff : this->cutoffs) {
          std::tuple<double> a1{cutoff};
          this->factory_args.emplace_back(filename, a1);
        }
      }
    }

    ~MultipleStructureManagerNLFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};
    const std::vector<double> cutoffs{{1., 2., 3.}};

    std::vector<Factory_t> factory_args{};
  };

  struct MultipleStructureManagerNLStrictFixture {
    using Factory_t = std::tuple<std::string,std::tuple<double>,std::tuple<double>>;
    using AdaptorTypeHolder_t = AdaptorTypeHolder<StructureManagerCenters, AdaptorNeighbourList, AdaptorStrict>;

    MultipleStructureManagerNLStrictFixture() {
      for (auto&& filename : this->filenames) {
        for (auto&& cutoff : this->cutoffs) {
          std::tuple<double> a1{cutoff};
          this->factory_args.emplace_back(filename, a1, a1);
        }
      }
    }

    ~MultipleStructureManagerNLStrictFixture() = default;

    const std::vector<std::string> filenames{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
        "reference_data/simple_cubic_8.json",
        "reference_data/small_molecule.json"};
    const std::vector<double> cutoffs{{1., 2., 3.}};

    std::vector<Factory_t> factory_args{};
  };

  template <class BaseFixture, class StructureManager, template<class> class ... AdaptorImplementationPack>
  struct MultipleStructureFixture : BaseFixture {
    using Parent = BaseFixture;
    using Factory_t = typename Parent::Factory_t;
    using Manager_t = typename internal::AdaptorTypeStacker<StructureManager,AdaptorImplementationPack...>::type;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    MultipleStructureFixture() : Parent{} {
      for (auto factory_arg : this->factory_args) {
        auto manager{call_with_tuple<Factory_t>::make_manager_stack(factory_arg)};
        this->managers.push_back(manager);
      }
    }

    template<typename T>
    struct call_with_tuple;

    template<typename ...T>
    struct call_with_tuple<std::tuple<T...>> {
      static ManagerPtr_t make_manager_stack(std::tuple<T...> tuple) {
        return helper(tuple, std::index_sequence_for<T...>());
      }
     protected:
      template<std::size_t... Is>
      static ManagerPtr_t helper(std::tuple<T...> tuple, std::index_sequence<Is...>) {
        ManagerPtr_t manager{make_structure_manager_stack<StructureManager, AdaptorImplementationPack...>(std::get<Is>(tuple)...)};
        return manager;
      }

    };


    ~MultipleStructureFixture() = default;

    std::list<std::shared_ptr<Manager_t>> managers{};
  };

}  // namespace rascal

#endif /* TEST_NEIGHBOURHOOD_H */
