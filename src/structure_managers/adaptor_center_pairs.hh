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

#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/cluster_ref_key.hh"

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
  struct StructureManager_traits<AdaptorCenterPairs<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
        ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    // New MaxOrder upon construction stays the same
    // todo(markus) reset PairOrder to zero since new pairs are added
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};

    using LayerByOrder = typename LayerExtender<
        MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
  };

  /**
   * Adaptor that adds pairs of centers, i.e. `ii` pairs to the list of
   * neighbours.
   */
  template <class ManagerImplementation>
  class AdaptorCenterPairs
      : public AdaptorFilter<ManagerImplementation,
                             ManagerImplementation::traits::MaxOrder> {
   public:
    using Manager_t = AdaptorCenterPairs<ManagerImplementation>;
    // using Parent = StructureManager<Manager_t>;
    using Parent = AdaptorFilter<ManagerImplementation,
                                 ManagerImplementation::traits::MaxOrder>;
    using traits = StructureManager_traits<AdaptorCenterPairs>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using Hypers_t = typename Parent::Hypers_t;

    //! Default constructor
    AdaptorCenterPairs() = delete;

    AdaptorCenterPairs(ImplementationPtr_t manager) : Parent{manager} {}

    AdaptorCenterPairs(ImplementationPtr_t manager, std::tuple<> /*val*/)
        : AdaptorCenterPairs{manager} {}

    AdaptorCenterPairs(ImplementationPtr_t manager,
                       const Hypers_t & /*adaptor_hypers*/)
        : AdaptorCenterPairs{manager} {}

    //! Copy constructor
    AdaptorCenterPairs(const AdaptorCenterPairs & other) = delete;

    //! Move constructor
    AdaptorCenterPairs(AdaptorCenterPairs && other) = delete;

    //! Destructor
    virtual ~AdaptorCenterPairs() = default;

    //! Copy assignment operator
    AdaptorCenterPairs & operator=(const AdaptorCenterPairs & other) = delete;

    //! Move assignment operator
    AdaptorCenterPairs & operator=(AdaptorCenterPairs && other) = delete;

    template <class... Args>
    void update(Args &&... arguments) {
      this->manager->update(std::forward<Args>(arguments)...);
    }

    //! This is where the magic happens and the center-pairs are added along
    //! with all previous pairs
    void perform_filtering() final {
      for (auto && atom : this->manager) {
        // get a ClusterRef<2> to construct an ii pair
        auto && ii_it{atom.begin()};
        // dereference for correct type
        auto && ii_pair{*ii_it};
        // reset atom tag to actually be an ii-pair
        ii_pair.set_atom_tag(1, atom.get_atom_tag());

        // add it to the list of pairs
        this->add_cluster(ii_pair);

        // add all other pairs in neighbour list
        for (auto && pair : atom) {
          this->add_cluster(pair);
        }
      }
    }
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_
