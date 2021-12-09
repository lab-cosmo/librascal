/**
 * @file test_manager_collection.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 June 2019
 *
 * @brief contains the fixtures with data for the tests of the structure
 *        manager collection
 *
 * @section LICENSE
 *
 * Copyright  2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef TESTS_TEST_MANAGER_COLLECTION_HH_
#define TESTS_TEST_MANAGER_COLLECTION_HH_

#include "test_adaptor.hh"

#include "rascal/structure_managers/structure_manager_collection.hh"

namespace rascal {

  struct StrictNLCCCollectionFixture
      : MultipleStructureManagerNLCCStrictFixture {
    using Parent = MultipleStructureManagerNLCCStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    StrictNLCCCollectionFixture() : Parent{} {};

    ~StrictNLCCCollectionFixture() = default;

    std::string filename{"reference_data/inputs/small_molecules-20.json"};
    int start{5};
    int length{3};
  };

  /**
   * BaseFixture is expected to be similar to
   * StrictNLCCCollectionFixture for instance
   */
  template <class BaseFixture>
  struct CollectionFixture : BaseFixture {
    using Parent = BaseFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
    using Manager_t = typename ManagerTypeHolder_t::type;
    using ManagerCollection_t =
        typename TypeHolderInjector<ManagerCollection, ManagerTypeList_t>::type;

    CollectionFixture() : Parent{} {
      for (auto & args : this->factory_args) {
        collections.emplace_back(args["adaptors"]);
        structures.emplace_back(args["structure"]);
      }
    }

    ~CollectionFixture() = default;

    json structures{};

    std::vector<ManagerCollection_t> collections{};
  };

}  // namespace rascal

#endif  // TESTS_TEST_MANAGER_COLLECTION_HH_
