/**
 * file test_manager_collection.hh
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

#ifndef TESTS_TEST_MANAGER_COLLECTION_HH_
#define TESTS_TEST_MANAGER_COLLECTION_HH_

#include "tests.hh"
#include "test_adaptor.hh"
#include "structure_managers/structure_manager_collection.hh"

namespace rascal {

  struct StrictNLCollectionFixture : MultipleStructureManagerNLCCStrictFixture {
    using Parent = MultipleStructureManagerNLCCStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    StrictNLCollectionFixture() : Parent{} {};

    ~StrictNLCollectionFixture() = default;

    std::string filename{"reference_data/dft-smiles_500.ubjson"};
    int start{5};
    int length{3};
  };

  /**
   * BaseFixture is expected to be similar to
   * StrictNLCollectionFixture for instance
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
