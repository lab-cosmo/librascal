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


namespace rascal {


  struct MultipleStructureSortedCoulomb
  {
    MultipleStructureSortedCoulomb() = default;
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<std::string> filenames{
      // "reference_data/CaCrP2O7_mvc-11955_symmetrized_.json",
      "reference_data/simple_cubic_8_.json"};
    std::vector<double> cutoffs{{1.9}};
    
    double central_decay{0.5};
    double interaction_cutoff{10};
    double interaction_decay{0.5};
    size_t size{15};
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