/**
 * file test_manager_collection.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   14 June 2019
 *
 * @brief
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

  struct ManagerNLCollectionFixture {
    using ManagerCollection_t = ManagerCollection<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    using Manager_t = typename ManagerCollection_t::Manager_t;

    ManagerNLCollectionFixture() {
      json ad1{
          {"name", "AdaptorNeighbourList"},
          {"initialization_arguments",
            {{"cutoff", cutoff},
            {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
      json ad2{{"name", "AdaptorStrict"},
                {"initialization_arguments", {{"cutoff", cutoff}}}};
      adaptors.emplace_back(ad1);
      adaptors.emplace_back(ad2);

      collection.set_adaptor_inputs(adaptors);
    }

    ~ManagerNLCollectionFixture() = default;
    const bool consider_ghost_neighbours{false};
    std::string filename{
        "reference_data/dft-smiles_500.ubjson"};
    const double cutoff{3.};

    json adaptors{};
    ManagerCollection_t collection{};

  };

}  // namespace rascal

#endif  // TESTS_TEST_MANAGER_COLLECTION_HH_
