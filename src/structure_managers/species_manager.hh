/**
 * file    species_manager.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   14 Sep 2018
 *
 * @brief Manager for expanding a structure manager into a collection
 *        of structure managers with separate combinations of species,
 *        class comment for SpeciesManager
 *
 * Copyright © 2018 Markus Stricker, Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_SPECIES_MANAGER_HH_
#define SRC_STRUCTURE_MANAGERS_SPECIES_MANAGER_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/property.hh"

#include <type_traits>
#include <array>
#include <tuple>
#include <map>
#include <memory>

namespace rascal {

  namespace internal {

    namespace detail {
      template <size_t Order>
      using Key_t = std::array<int, Order>;
      template <size_t Order, class SpeciesManager>
      using Value_t =
        std::unique_ptr<typename SpeciesManager::template Filter<Order>>;
      template <size_t Order, class SpeciesManager>
      using Map_t =
          std::map<Key_t<Order>, Value_t<Order, SpeciesManager>>;
    }  // namespace detail

    template <class SpeciesManager, size_t... OrdersMinusOne>
    auto get_filter_container_helper(
        std::index_sequence<OrdersMinusOne...> /*orders*/) -> decltype(auto) {
      return std::tuple<
          detail::Map_t<OrdersMinusOne + 1, SpeciesManager>...>{};
    }

    template <class SpeciesManager, size_t MaxOrder>
    auto get_filter_container() -> decltype(auto) {
      return get_filter_container_helper<SpeciesManager>(
          std::make_index_sequence<MaxOrder>{});
    }

    template <class SpeciesManager, size_t MaxOrder>
    using FilterContainer_t =
        decltype(get_filter_container<SpeciesManager, MaxOrder>());

  }  // namespace internal

  /**
   * Takes a Structure manager and splits it into sub sets
   * distinguished by species combinations see illustration for a
   * example with 2 species, a and b:
   *
   * ```
   *               SpeciesManager
   *               /             \
   *              a               b
   *            /    \          /    \
   *          aa      ab      ba      bb
   *         /   \   /   \   /   \   /   \
   *        aaa aab aba abb baa bab bba bbb
   * ```
   *
   * Use case example, we have a function `fun` to evaluate on all triplets
   * of type aba:
   *
   * ```
   * SpeciesManager species_manager{...};
   * std::array<int, 3> species_indices{a, b, a};
   * fun(species_manager[species_indices])
   * ```
   */
  template <class ManagerImplementation, size_t MaxOrder>
  class SpeciesManager {
   public:
    using traits = StructureManager_traits<ManagerImplementation>;
    template <size_t NbElements>
    using SpeciesMap_t =
        internal::detail::Map_t<NbElements, SpeciesManager>;

    /**
     * implementation of AdaptorFilter for species filtering
     */
    template <size_t Order>
    class Filter;

    using FilterContainer_t =
        internal::FilterContainer_t<SpeciesManager, MaxOrder>;

    static_assert(traits::MaxOrder <= MaxOrder,
                  "MaxOrder of underlying manager is insufficient.");

    //! Default constructor
    SpeciesManager() = delete;

    explicit SpeciesManager(ManagerImplementation & manager);

    //! Copy constructor
    SpeciesManager(const SpeciesManager & other) = delete;

    //! Move constructor
    SpeciesManager(SpeciesManager && other) = default;

    //! Destructor
    virtual ~SpeciesManager() = default;

    //! Copy assignment operator
    SpeciesManager & operator=(const SpeciesManager & other) = delete;

    //! Move assignment operator
    SpeciesManager & operator=(SpeciesManager && other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes building either the neighbour list or to
     * make triplets, quadruplets, etc. depending on the MaxOrder
     */
    void update();

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments) {
      this->update();
      this->structure_manager.update(arguments...);
    }

    template <size_t Order>
    Filter<Order> &
    operator[](const std::array<int, Order> & species_indices) {
      auto & filter_map{std::get<Order - 1>(this->filters)};
      auto location{filter_map.find(species_indices)};

      if (location == filter_map.end()) {
        // this species combo is not yet in the container, therefore
        // create new empty one
        auto new_filter{
            std::make_unique<Filter<Order>>(this->structure_manager)};
        // insertion returns a ridiculous type, see spec
        auto new_location{std::get<0>(
            filter_map.emplace(species_indices, std::move(new_filter)))};
        return *std::get<1>(*new_location);
      } else {
        return *std::get<1>(*location);
      }
    }

    template <size_t NbElements>
    SpeciesMap_t<NbElements> & filters_by_nb_elements() {
      return std::get<NbElements>(this->filters);
    }

   protected:
    //! underlying structure manager to be filtered upon update()
    ManagerImplementation & structure_manager;
    //! storage by cluster order for the filtered managers
    FilterContainer_t filters;
  };

  template <class ManagerImplementation, size_t MaxOrder>
  template <size_t Order>
  class SpeciesManager<ManagerImplementation, MaxOrder>::Filter
      : public AdaptorFilter<ManagerImplementation, Order>
  {
   public:

    using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;
    using Parent = AdaptorFilter<ManagerImplementation, Order>;
    //! Default constructor
    Filter() = delete;

    Filter(SpeciesManager_t & species_manager)
        : Parent{species_manager.get_manager()},
          species_manager{species_manager} {}

    //! Copy constructor
    Filter(const Filter &other) = delete;

    //! Move constructor
    Filter(Filter &&other) = delete;

    //! Destructor
    virtual ~Filter() = default;

    //! Copy assignment operator
    Filter& operator=(const Filter &other) = delete;

    //! Move assignment operator
    Filter& operator=(Filter &&other) = delete;

    void perform_filtering() final {};

   protected:
    SpeciesManager_t& species_manager;
  };

  /* ----------------------------------------------------------------------
       */
      template <class ManagerImplementation, size_t MaxOrder>
      SpeciesManager<ManagerImplementation, MaxOrder>::SpeciesManager(
          ManagerImplementation & manager)
      : structure_manager{manager},
        filters{
            internal::get_filter_container<ManagerImplementation, MaxOrder>()} {
  }

  namespace internal {
    /**
     * Helper struct that loops over a cluster or manager, and
     * segregates the iteratee (i.e. the next higher order clusters)
     * by species. If the loop has not reached the highest cluster
     * order (i.e. MaxOrder), it recursively loops also over the next
     * higher order clusters.
     */
    template <class ManagerImplementation, size_t MaxOrder,
              size_t Remaining = MaxOrder>
    struct FilterSpeciesLoop {
      using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;
      using NextFilterSpeciesLoop =
          FilterSpeciesLoop<ManagerImplementation, MaxOrder, Remaining - 1>;

      template <class Cluster>
      static void loop(Cluster & cluster, SpeciesManager_t & species_manager) {
        constexpr auto ClusterOrder{MaxOrder-Remaining + 1};
        // reset all filters
        for (auto && tup :
             species_manager.template filters_by_nb_elements<ClusterOrder>()) {
          static_assert(tup.first.size() == ClusterOrder,
                        "FilterSpeciesLoop constructed with wrong template "
                        "parameters");
          tup.second.reset_initial_state();
        }
        // refill all filters
        for (auto && next_cluster : cluster) {
          auto && species_indices{next_cluster.get_atom_types()};
          species_manager[species_indices].add_cluster(next_cluster);
          NextFilterSpeciesLoop::loop(next_cluster, species_manager);
        }
      }
    };

    /**
     * Recursion tail of the helper loop which does nothing at all
     */
    template <class ManagerImplementation, size_t MaxOrder>
    struct FilterSpeciesLoop<ManagerImplementation, MaxOrder, 0> {
      using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;
      template <class Cluster>
      static void loop(Cluster & /*cluster*/,
                       SpeciesManager_t & /*species_manager*/) {}
    };

  }  // namespace internal

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation, size_t MaxOrder>
  void SpeciesManager<ManagerImplementation, MaxOrder>::update() {
    // number of levels to be descended into is know at compile time,
    // but not at writing time, hence the indirection to
    // FilterSpeciesLoop
    using FilterSpeciesLoop =
        internal::FilterSpeciesLoop<ManagerImplementation, MaxOrder>;
    FilterSpeciesLoop::loop(this->structure_manager, *this);
  }

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_SPECIES_MANAGER_HH_
