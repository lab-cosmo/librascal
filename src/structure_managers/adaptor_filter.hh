/**
 * file   adaptor_filter.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   23 Oct 2018
 *
 * @brief An adaptor that provides a filtered (masked) view
 *        on an existing structure manager.
 *
 * Copyright © 2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef ADAPTOR_FILTER_H
#define ADAPTOR_FILTER_H

#include <type_traits>

namespace rascal {

  /**
   * Forward declaration for traits
   */
  template <class StructureManagerImplementation, size_t MaxOrder>
  class AdaptorFilter;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation, size_t MaxOrder>
  struct StructureManager_traits
  <AdaptorFilter<ManagerImplementation, MaxOrder>> {
    using traits =
      StructureManager_traits<AdaptorImplementation::ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{traits::Strict};
    constexpr static bool HasDistances{traits::HasDistances};
    constexpr static bool HasDirectionVectors{traits::HasDirectionVectors};
    constexpr static int Dim{traits::Dim};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{MaxOrder};
    //! New Layer
    //! TODO: Is this the correct way to initialize the increased order?
    using LayerByOrder =
      typename LayerIncreaser
      <MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Adaptor that increases the MaxOrder of an existing StructureManager. This
   * means, if the manager does not have a neighbourlist, there is nothing this
   * adaptor can do (hint: use adaptor_neighbour_list before and stack this on
   * top), if it exists, triplets, quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation, size_t MaxOrder>
  class AdaptorFilter: public
  StructureManager<AdaptorFilter<ManagerImplementation>> {
  public:
    using Parent = StructureManager<AdaptorFilter<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorFilter>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;

    static_assert(traits::MaxOrder <= ManagerImplementation::traits::MaxOrder,
                  "can only present view on existing clusters");

    //! Default constructor
    AdaptorFilter() = delete;

    //! constructor underlying manager
    explicit AdaptorFilter(ManagerImplementation & manager): manager{manager} {}

    //! Copy constructor
    AdaptorFilter(const AdaptorFilter &other) = delete;

    //! Move constructor
    AdaptorFilter(AdaptorFilter &&other) = default;

    //! Destructor
    virtual ~AdaptorFilter() = default;

    //! Copy assignment operator
    AdaptorFilter& operator=(const AdaptorFilter &other) = delete;

    //! Move assignment operator
    AdaptorFilter& operator=(AdaptorFilter &&other) = default;

    //! returns the distance between atoms in a given pair
    template <size_t Order, size_t Layer>
    inline const std::enable_if_t<traits::HasDistances, double> &
    get_distance(const ClusterRefKey<Order, Layer> &
                                       pair) const {
      return this->distance[pair];
    }


  protected:
    ManagerImplementation & manager;
  private:
  };





}  // rascal

#endif /* ADAPTOR_FILTER_H */
