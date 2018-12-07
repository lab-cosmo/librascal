/**
 * file   test_representation_manager_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TEST_REPRESENTATION_H
#define TEST_REPRESENTATION_H

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"
#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "json_io.hh"
#include "rascal_utility.hh"

#include <tuple>

namespace rascal {


  struct MultipleStructureSortedCoulomb {
    MultipleStructureSortedCoulomb() {}
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<std::string> filenames{
      "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
      "reference_data/simple_cubic_8.json",
      "reference_data/small_molecule.json"
      };
    std::vector<double> cutoffs{{1., 2., 3.}};

    std::list<json> hypers{
      {{"central_decay", 0.5},
      {"interaction_cutoff", 10.},
      {"interaction_decay", 0.5},
      {"size", 120}}
      };
  };

  struct SortedCoulombTestData {
    SortedCoulombTestData() {
      std::vector<std::uint8_t> ref_data_ubjson;
      internal::read_binary_file(this->ref_filename, ref_data_ubjson);
      json ref_data = json::from_ubjson(ref_data_ubjson);
      filenames = ref_data.at("filenames").get<std::vector<std::string>>();
      cutoffs = ref_data.at("cutoffs").get<std::vector<double>>();
      data_sort_distance = ref_data.at("distance").get<json>();
      data_sort_rownorm = ref_data.at("rownorm").get<json>();
    }
    ~SortedCoulombTestData() = default;

    std::string ref_filename{"reference_data/sorted_coulomb_reference.ubjson" };
    std::vector<std::string> filenames{};
    std::vector<double> cutoffs{};
    json data_sort_distance{};
    json data_sort_rownorm{};
    json hypers{};
    json feature_matrices{};
  };

  template< class StructureManager,
            template<typename, Option ...opts > class RepresentationManager,
            class BaseFixture, Option ...options_>
  struct RepresentationFixture
  :MultipleStructureManagerStrictFixture<StructureManager, BaseFixture> {
    using Parent = MultipleStructureManagerStrictFixture<StructureManager,
                                                         BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = RepresentationManager<Manager_t, options_...>;

    RepresentationFixture() = default;
    ~RepresentationFixture() = default;

    std::list<Representation_t> representations{};
    std::vector<Option> options{options_...};
  };

/* ---------------------------------------------------------------------- */

} // RASCAL

#endif /* TEST_REPRESENTATION_H */
