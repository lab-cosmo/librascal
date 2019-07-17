/**
 * file   adaptor_center_pairs.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *         Markus Stricker <markus.stricker@epfl.ch>
 * @date   12 July 2019
 *
 * @brief implements an adaptor for structure_managers, adding the central atom
 * to the iteration over the neighbours
 *
 * Copyright  2019 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/updateable_base.hh"

namespace rascal {
  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorCenterPairs;

  /**
   * Specialisation of traits for reset pair <code>MaxOrder</code> for adding
   * center-center pairs
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorMaxOrder<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
        ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    // New MaxOrder upon construction stays the same
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};

    using LayerByOrder = typename LayerExtender<
        MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
  };

}  // namespace rascal

#endif // SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_
