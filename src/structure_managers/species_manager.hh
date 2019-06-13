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
 * Copyright  2018 Markus Stricker, Till Junge, COSMO (EPFL), LAMMM (EPFL)
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
#include "structure_managers/updateable_base.hh"
#include "utils/key_standardisation.hh"

#include <type_traits>
#include <array>
#include <tuple>
#include <map>
#include <memory>

namespace rascal {

  namespace internal {

    namespace detail {
      template <size_t Order>
      using Key_t = KeyStandardisation<int, Order>;
      using Value_t = std::unique_ptr<FilterBase>;
      template <size_t Order>
      using Map_t = std::map<Key_t<Order>, Value_t>;
    }  // namespace detail
  }    // namespace internal

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
  class SpeciesManager : public Updateable,
                         public std::enable_shared_from_this<
                             SpeciesManager<ManagerImplementation, MaxOrder>> {
   public:
    using traits = StructureManager_traits<ManagerImplementation>;
    using Manager_t = SpeciesManager<ManagerImplementation, MaxOrder>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using Key_t = internal::detail::Key_t<MaxOrder>;

    /**
     * implementation of AdaptorFilter for species filtering
     */
    template <size_t Order>
    class Filter;

    using FilterContainer_t = typename internal::detail::Map_t<MaxOrder>;

    static_assert(traits::MaxOrder <= MaxOrder,
                  "MaxOrder of underlying manager is insufficient.");

    //! Default constructor
    SpeciesManager() = delete;

    explicit SpeciesManager(ImplementationPtr_t manager);

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
    void update_self();

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments) {
      // SpeciesManager does not belong to the chain of updatables
      // like a regular adaptor so you need to set the status here
      if (sizeof...(arguments) > 0) {
        this->set_update_status(false);
      }
      this->structure_manager->update(std::forward<Args>(arguments)...);
      /**
       * Note that the update and update_self rational is hacked here for
       * efficiency.
       * The species_manager::Filter don't know how to update themselves
       * because this functionality is defered to species_manager. For example
       * accros a MD simulation the possible atom types might change and
       * SpeciesManager can change dynamically the list of filters but they
       * can't do it themselves.
       */
      this->update_self();
    }

    /**
     * When the underlying structure changes, all computations are potentially
     * invalid. This function triggers the setting of the statue variable to
     * `false` along the tree to the managers and the properties it holds.
     */
    void send_changed_structure_signal() final {
      this->set_update_status(false);
      for (auto && child : this->children) {
        if (not child.expired()) {
          child.lock()->send_changed_structure_signal();
        }
      }
    }

    /**
     * This function is never called because SpeciesManager belongs to the
     * update upward tree but only the this->filters belong to the downward
     * tree that would typically trigger the update_self.
     */
    void update_children() final {}

    template <size_t Order>
    Filter<Order> & operator[](const std::array<int, Order> & species_indices) {
      auto && location{this->filters.find(Key_t{species_indices})};

      if (location == this->filters.end()) {
        // this species combo is not yet in the container, therefore
        // create new empty one
        auto new_filter{std::make_unique<Filter<Order>>(*this)};
        // insertion returns a ridiculous type: ((Key, Value), success), where
        // Value is the filter
        // create a (species_indices, filter) pair and add it to the list
        auto && retval{this->filters.emplace(Key_t{species_indices},
                                             std::move(new_filter))};
        auto && new_location{retval.first};
        return static_cast<Filter<Order> &>(*(new_location->second));
      } else {
        return static_cast<Filter<Order> &>(*(location->second));
      }
    }

    ImplementationPtr_t get_structure_manager() const {
      return this->structure_manager;
    }

   protected:
    //! underlying structure manager to be filtered upon update_self()
    ImplementationPtr_t structure_manager;
    //! storage by cluster order for the filtered managers
    FilterContainer_t filters{};
  };

  template <class ManagerImplementation, size_t MaxOrder>
  template <size_t Order>
  class SpeciesManager<ManagerImplementation, MaxOrder>::Filter
      : public AdaptorFilter<ManagerImplementation, Order> {
   public:
    using SpeciesManager_t = SpeciesManager<ManagerImplementation, MaxOrder>;
    using Parent = AdaptorFilter<ManagerImplementation, Order>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;

    //! Default constructor
    Filter() = delete;

    Filter(SpeciesManager_t & species_manager)
        : Parent{species_manager.get_structure_manager()},
          species_manager{species_manager} {}

    //! Copy constructor
    Filter(const Filter & other) = delete;

    //! Move constructor
    Filter(Filter && other) = delete;

    //! Destructor
    virtual ~Filter() = default;

    //! Copy assignment operator
    Filter & operator=(const Filter & other) = delete;

    //! Move assignment operator
    Filter & operator=(Filter && other) = delete;

    void perform_filtering() final{};

   protected:
    SpeciesManager_t & species_manager;
  };

  /* ----------------------------------------------------------------------
   */
  template <class ManagerImplementation, size_t MaxOrder>
  SpeciesManager<ManagerImplementation, MaxOrder>::SpeciesManager(
      ImplementationPtr_t manager)
      : structure_manager{manager} {}

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
  void SpeciesManager<ManagerImplementation, MaxOrder>::update_self() {
    // reset all filters before filling them again
    for (auto && key_filter : this->filters) {
      key_filter.second->reset_initial_state();
    }

    // number of levels to be descended into is know at compile time,
    // but not at writing time, hence the indirection to
    // FilterSpeciesLoop
    using FilterSpeciesLoop =
        internal::FilterSpeciesLoop<ManagerImplementation, MaxOrder>;
    FilterSpeciesLoop::loop(this->structure_manager, *this);
  }

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_SPECIES_MANAGER_HH_
