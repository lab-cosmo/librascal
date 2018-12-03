/**
 * file   test_representation_manager_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright © 2018 Musil Felix, Max Veit, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
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
#include "representations/representation_manager_spherical_expansion.hh"


namespace rascal {


  struct MultipleStructureSortedCoulomb
  {
    MultipleStructureSortedCoulomb() = default;
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<std::string> filenames{
      "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
      "reference_data/simple_cubic_8.json",
      "reference_data/small_molecule.json"
      };
    std::vector<double> cutoffs{{1,2,3}};

    std::list<json> hypers{
      {{"central_decay",0.5},
      {"interaction_cutoff",10},
      {"interaction_decay",0.5},
      {"size",120}}
      };
  };

  struct MultipleStructureSphericalExpansion
  {
    MultipleStructureSphericalExpansion() = default;
    ~MultipleStructureSphericalExpansion() = default;

    std::vector<std::string> filenames {
      "reference_data/simple_cubic_8.json",
      "reference_data/small_molecule.json"
    };
    std::vector<double> cutoffs{{1,2,3}};

    std::list<json> hypers{ {
      {"interaction_cutoff", 6.0},
      {"cutoff_smooth_width", 1.0},
      {"max_radial", 10},
      {"max_angular", 8},
      {"gaussian_sigma_type", "Constant"},
      {"gaussian_sigma", 0.5}
    }};
  };

  template< class StructureManager,
            template<typename> class RepresentationManager,
            class BaseFixture>
  struct RepresentationFixture
  :MultipleStructureManagerStrictFixture<StructureManager, BaseFixture>
  {
    using Parent = MultipleStructureManagerStrictFixture<StructureManager,
                    BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = RepresentationManager<Manager_t>;

    RepresentationFixture() = default;
    ~RepresentationFixture() = default;

    std::list<Representation_t> representations{};
  };

/* ---------------------------------------------------------------------- */

} // RASCAL

#endif /* TEST_REPRESENTATION_H */
