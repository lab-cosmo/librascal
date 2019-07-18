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
    //using Parent = StructureManager<Manager_t>;
    using Parent = AdaptorFilter<ManagerImplementation,
                                 ManagerImplementation::traits::MaxOrder>;
    using traits = StructureManager_traits<AdaptorCenterPairs>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using Hypers_t = typename Parent::Hypers_t;

    //! Default constructor
    AdaptorCenterPairs() = delete;

    AdaptorCenterPairs(ImplementationPtr_t manager)
      : Parent{manager} {}

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

    //! This is where the magic happens and the center-pairs are added
    void perform_filtering() final {
      for (auto && atom : this->manager) {

        // construct and add ii-pair
        // auto ii_pair{*(atom.begin())};

        // ii_pair.back() = ii_pair.front();
        // ii_pair.get_atom_type() = atom.get_atom_type();
        // ii_pair.get_position() = atom.get_position();

        // this->add_cluster(ii_pair);

        // add all other pairs in neighbour list
        for (auto && pair : atom) {
          this->add_cluster(pair);
        }
      }
    }
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_

//   StructureManager<AdaptorCenterPairs<ManagerImplementation>>,
//         public std::enable_shared_from_this<
//             AdaptorCenterPairs<ManagerImplementation>> {

//    public:
//     using Manager_t = AdaptorCenterPairs<ManagerImplementation>;
//     using Parent = StructureManager<Manager_t>;
//     using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
//     using traits = StructureManager_traits<AdaptorCenterPairs>;
//     using Vector_ref = typename Parent::Vector_ref;
//     using Hypers_t = typename Parent::Hypers_t;

//     static_assert(traits::MaxOrder > 1,
//                   "ManagerImplementation needs to handle pairs");
//     constexpr static auto AtomLayer{
//         Manager_t::template cluster_layer_from_order<1>()};
//     constexpr static auto PairLayer{
//         Manager_t::template cluster_layer_from_order<2>()};

//     class Filter;

//     //! Default constructor
//     AdaptorCenterPairs() = delete;

//     //! Construct a neighbourhood which includes pairs of centers
//     AdaptorCenterPairs(ImplementationPtr_t manager)
//         : structure_manager{manager} {};

//     //! Copy constructor
//     AdaptorCenterPairs(const AdaptorCenterPairs & other) = delete;

//     //! Move constructor
//     AdaptorCenterPairs(AdaptorCenterPairs && other) = default;

//     //! Destructor
//     virtual ~AdaptorCenterPairs() = default;

//     //! Copy assignment operator
//     AdaptorCenterPairs & operator=(const AdaptorCenterPairs & other) =
//     delete;

//     //! Move assignment operator
//     AdaptorCenterPairs & operator=(AdaptorCenterPairs && other) = default;

//     //! Update the adaptor assuming the underlying manager was update
//     inline void update_self();

//     //! Update the underlying manager as well as the adaptor
//     template <class... Args>
//     void update(Args &&... arguments);

//    protected:
//     ImplementationPtr_t structure_manager;

//    private:
//   };

//   template <class ManagerImplementation>
//   class AdaptorCenterPairs<ManagerImplementation>::Filter
//       : public AdaptorFilter<ManagerImplementation, traits::MaxOrder> {
//    public:
//     using Parent = AdaptorFilter<ManagerImplementation, traits::MaxOrder>;
//     using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;

//     //! Default constructor
//     Filter() = delete;

//     Filter(ImplementationPtr_t manager)
//         : Parent{manager}, structure_manager{manager} {}

//     //! Copy constructor
//     Filter(const Filter & other) = delete;

//     //! Move constructor
//     Filter(Filter && other) = delete;

//     //! Destructor
//     virtual ~Filter() = default;

//     //! Copy assignment operator
//     Filter & operator=(const Filter & other) = delete;

//     //! Move assignment operator
//     Filter & operator=(Filter && other) = delete;

//     //! This is where the magic happens and the center-pairs are added
//     void perform_filtering() final {
//       for (auto && atom : this->structure_manager) {

//         // construct and add ii-pair
//         auto ii_pair{*(atom.begin())};

//         ii_pair.back() = ii_pair.front();
//         ii_pair.get_atom_type() = atom.get_atom_type();
//         ii_pair.get_position() = atom.get_position();

//         this->add_cluster(ii_pair);

//         // add all other pairs in neighbour list
//         for (auto && pair : atom) {
//           this->add_cluster(pair);
//         }
//       }
//     };

//    protected:
//     ImplementationPtr_t structure_manager;
//   };
