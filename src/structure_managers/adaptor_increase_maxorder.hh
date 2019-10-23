/**
 * @file   adaptor_increase_maxorder.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jun 2018
 *
 * @brief implements an adaptor for structure_managers, which
 * creates a full and half neighbourlist if there is none and
 * triplets/quadruplets, etc. if existent.
 *
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_

#include "basic_types.hh"
#include "lattice.hh"
#include "rascal_utility.hh"
#include "structure_managers/property.hh"
#include "structure_managers/structure_manager.hh"

#include <set>
#include <vector>
#include <algorithm>

namespace rascal {
  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorMaxOrder<ManagerImplementation>> {
    using parent_traits = StructureManager_traits<ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
        parent_traits::HasDirectionVectors};
    constexpr static int Dim{parent_traits::Dim};
    constexpr static bool HasCenterPair{parent_traits::HasCenterPair};
    constexpr static int StackLevel{parent_traits::StackLevel + 1};
    // New MaxOrder upon construction
    constexpr static size_t MaxOrder{parent_traits::MaxOrder + 1};
    // Extend the layer by one with the new MaxOrder
    using LayerByOrder =
        typename LayerExtender<MaxOrder,
                               typename parent_traits::LayerByOrder>::type;
    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
        parent_traits::NeighbourListType};
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Adaptor that increases the MaxOrder of an existing StructureManager. This
   * means, if the manager does not have a neighbourlist, there is nothing this
   * adaptor can do (hint: use adaptor_neighbour_list before and stack this on
   * top), if it exists, triplets, quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder
      : public StructureManager<AdaptorMaxOrder<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorMaxOrder<ManagerImplementation>> {
   public:
    using Manager_t = AdaptorMaxOrder<ManagerImplementation>;
    using Parent = StructureManager<Manager_t>;
    using ManagerImplementation_t = ManagerImplementation;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorMaxOrder>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    static_assert(traits::MaxOrder > 2,
                  "ManagerImplementation needs at least a pair list for"
                  " extension.");

    //! Default constructor
    AdaptorMaxOrder() = delete;

    /**
     * Given at least a pair list, this adaptor creates the next Order
     * list. I.e. from pairs to triplets, triplets to quadruplet, etc. Not
     * cutoff is needed, the cutoff is implicitly given by the neighbourlist,
     * which was built
     */
    explicit AdaptorMaxOrder(ImplementationPtr_t manager);

    AdaptorMaxOrder(ImplementationPtr_t manager, std::tuple<>)
        : AdaptorMaxOrder(manager) {}

    AdaptorMaxOrder(ImplementationPtr_t manager,
                    const Hypers_t & /*adaptor_hypers*/)
        : AdaptorMaxOrder(manager) {}

    //! Copy constructor
    AdaptorMaxOrder(const AdaptorMaxOrder & other) = delete;

    //! Move constructor
    AdaptorMaxOrder(AdaptorMaxOrder && other) = default;

    //! Destructor
    virtual ~AdaptorMaxOrder() = default;

    //! Copy assignment operator
    AdaptorMaxOrder & operator=(const AdaptorMaxOrder & other) = delete;

    //! Move assignment operator
    AdaptorMaxOrder & operator=(AdaptorMaxOrder && other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes making triplets, quadruplets,
     * etc. depending on the MaxOrder, pair list has to be present.
     */
    void update_self();

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    bool get_consider_ghost_neighbours() const { return true; }

    /**
     * Returns the linear indices of the clusters (whose atom tags are stored
     * in counters). For example when counters is just the list of atoms, it
     * returns the index of each atom. If counters is a list of pairs of indices
     * (i.e. specifying pairs), for each pair of indices i,j it returns the
     * number entries in the list of pairs before i,j appears.
     */
    template <size_t Order>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const;

    //! Returns the number of clusters of size cluster_size
    size_t get_nb_clusters(size_t order) const {
      switch (order) {
      case traits::MaxOrder:
        return this->neighbours_atom_tag.size();
      default:
        return this->manager->get_nb_clusters(order);
      }
    }

    size_t get_size_with_ghosts() const {
      return this->manager->get_size_with_ghosts();
    }

    //! Returns number of clusters of the original manager
    size_t get_size() const { return this->manager->get_size(); }

    //! Returns position of an atom with index atom_tag
    Vector_ref get_position(size_t atom_tag) {
      return this->manager->get_position(atom_tag);
    }

    //! Returns position of the given atom object (useful for users)
    Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager->get_position(atom.get_index());
    }

    //! get atom type from underlying manager
    int get_atom_type(int atom_tag) const {
      return this->manager->get_atom_type(atom_tag);
    }

    //! return atom type
    int get_atom_type(const AtomRef_t & atom) const {
      return this->manager->get_atom_type(atom.get_atom_tag());
    }

    /**
     * Returns the id of the index-th (neighbour) atom of the cluster that is
     * the full structure/atoms object, i.e. simply the id of the index-th atom
     */
    int get_neighbour_atom_tag(const Parent &, size_t index) const {
      return this->manager->get_neighbour_atom_tag(*this->manager, index);
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template <size_t Order, size_t Layer>
    int get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & cluster,
                               size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");

      // necessary helper construct for static branching
      using IncreaseHelper_t =
          internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

      if (Order < (traits::MaxOrder - 1)) {
        return IncreaseHelper_t::get_neighbour_atom_tag(*this->manager, cluster,
                                                        index);
      } else {
        auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
        return this->neighbours_atom_tag[offset + index];
      }
    }

    size_t get_atom_index(const int atom_tag) const {
      return this->manager->get_atom_index(atom_tag);
    }

    //! Returns the number of neighbors of a given cluster
    template <size_t Order, size_t Layer>
    size_t
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");
      /*
       * Here it is traits::MaxOrder-1, because only the current manager has the
       * right answer to the number of neighbours of the MaxOrder-1 tuple. This
       * is the 'else' case.
       */

      // necessary helper construct for static branching
      using IncreaseHelper_t =
          internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

      if (Order < (traits::MaxOrder - 1)) {
        return IncreaseHelper_t::get_cluster_size(*this->manager, cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
        return this->nb_neigh[access_index];
      }
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager() {
      return this->manager->get_shared_ptr();
    }

   protected:
    template <bool IsCompactCluster>
    void update_self_helper();

    //! Extends the list containing the number of neighbours with a 0
    void add_entry_number_of_neighbours() { this->nb_neigh.push_back(0); }

    //! Adds a given atom tag as new cluster neighbour
    void add_neighbour_of_cluster(const int atom_tag) {
      // adds `atom_tag` to neighbours
      this->neighbours_atom_tag.push_back(atom_tag);
      // increases the number of neighbours
      this->nb_neigh.back()++;
    }

    //! Sets the correct offsets for accessing neighbours
    void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      if (n_tuples > 0) {
        this->offsets.reserve(n_tuples);
        this->offsets.resize(1);
        for (size_t i{0}; i < n_tuples; ++i) {
          this->offsets.emplace_back(this->offsets[i] + this->nb_neigh[i]);
        }
      }
    }

    //! reference to underlying manager
    ImplementationPtr_t manager;

    //! Construct for reaching the MaxOrder and adding neighbours of at MaxOrder
    //! (the work is done at the Order-th recursion, when IsTail is true)
    template <size_t Order, AdaptorTraits::NeighbourListType NeighbourListType,
              bool IsCompactCluster, bool IsTail>
    struct AddOrderLoop;

    //! Stores the number of neighbours for every traits::MaxOrder-1-clusters
    std::vector<size_t> nb_neigh{};

    //! Stores all neighbours atom tag of traits::MaxOrder-1-clusters
    std::vector<int> neighbours_atom_tag{};

    /**
     * Stores the offsets of traits::MaxOrder-1-*clusters for accessing
     * `neighbours`, from where nb_neigh can be counted
     */
    std::vector<size_t> offsets{};

    bool compute_compact_clusters{false};

   private:
  };

  /* ---------------------------------------------------------------------- */
  //! Constructor of the next level manager
  template <class ManagerImplementation>
  AdaptorMaxOrder<ManagerImplementation>::AdaptorMaxOrder(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, nb_neigh{},
        neighbours_atom_tag{}, offsets{} {
    if (traits::MaxOrder < 3) {
      throw std::runtime_error("Increase MaxOrder: No pair list in underlying"
                               " manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  //! update, involing the update of the underlying manager
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorMaxOrder<ManagerImplementation>::update(Args &&... arguments) {
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  //! structure for static looping up until pair order
  template <class ManagerImplementation>
  template <size_t Order, AdaptorTraits::NeighbourListType NeighbourListType,
            bool IsCompactCluster, bool IsTail>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;

    using NextOrderLoop =
        AddOrderLoop<Order + 1, NeighbourListType, IsCompactCluster,
                     (Order + 1 == OldMaxOrder)>;

    using AtomClusterRef_t =
        typename ManagerImplementation::template ClusterRef<1>;

    // do nothing, if MaxOrder is not reached, except call the next order
    static void loop(AtomClusterRef_t & atom, ClusterRef_t & cluster,
                     size_t start_index,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {
      size_t new_start_index{
          (start_index + 1) *
          static_cast<int>(NeighbourListType ==
                           AdaptorTraits::NeighbourListType::half)};

      using Iterator_t = typename ClusterRef_t::iterator;
      using NextCluster_t = typename Iterator_t::value_type;
      for (Iterator_t next_cluster_it{cluster.get_iterator_at(start_index)};
           next_cluster_it != cluster.end(); ++next_cluster_it) {
        auto && next_cluster{*next_cluster_it};
        constexpr size_t previous_order{NextCluster_t::order() - 1};
        auto & next_cluster_indices{
            std::get<previous_order>(manager.cluster_indices_container)};

        // keep copying underlying cluster indices, they are not changed
        auto indices{next_cluster.get_cluster_indices()};
        next_cluster_indices.push_back(indices);

        NextOrderLoop::loop(atom, next_cluster, new_start_index, manager);
        if (NeighbourListType == AdaptorTraits::NeighbourListType::half) {
          ++new_start_index;
        }
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * At desired MaxOrder (plus one), here is where the magic happens and the
   * neighbours of the same order are added as the Order+1.  add check for non
   * half neighbour list.
   *
   * TODO: currently, this implementation is not distinguishing between minimal
   * and full lists. E.g. this has to be adjusted to include both, the i- and
   * the j-atoms of each pair as an i-atom in a triplet (center).
   */
  template <class ManagerImplementation>
  template <size_t Order, AdaptorTraits::NeighbourListType NeighbourListType,
            bool IsCompactCluster>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop<
      Order, NeighbourListType, IsCompactCluster, true> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};

    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using AtomClusterRef_t =
        typename ManagerImplementation::template ClusterRef<1>;

    using traits = typename AdaptorMaxOrder<ManagerImplementation>::traits;

    //! loop through the orders to get to the maximum order, this is agnostic to
    //! the underlying MaxOrder, just goes to the maximum
    static void loop(AtomClusterRef_t & atom, ClusterRef_t & cluster,
                     size_t start_index,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {
      // add an entry for the current clusters' neighbours
      manager.add_entry_number_of_neighbours();

      // collect all possible neighbours of the cluster: collection of all
      // neighbours of current _central_ atoms for a centered cluster.
      if (NeighbourListType == AdaptorTraits::NeighbourListType::half) {
        for (auto pair_it{atom.get_iterator_at(start_index)};
             pair_it != atom.end(); ++pair_it) {
          if (IsCompactCluster) {
            throw std::runtime_error("Not implemented yet.");
          } else {
            manager.add_neighbour_of_cluster((*pair_it).back());
          }
        }
      } else {
        auto && i_atoms{cluster.get_atom_tag_list()};
        for (auto pair : atom) {
          auto && j_atom{pair.back()};
          if (std::find(i_atoms.begin(), i_atoms.end(), j_atom) ==
              i_atoms.end()) {
            if (IsCompactCluster) {
              throw std::runtime_error("Not implemented yet.");
            } else {
              manager.add_neighbour_of_cluster(j_atom);
            }
          }
        }
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Dispatch to the approprate loop
   */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::update_self() {
    if (this->compute_compact_clusters) {
      this->template update_self_helper<true>();
    } else {
      this->template update_self_helper<false>();
    }
  }
  /**
   * This is the loop, which runs recursively goes to the maximum Order and
   * then increases it by one (i.e. pairs->triplets, triplets->quadruplets,
   * etc.
   */
  template <class ManagerImplementation>
  template <bool IsCompactCluster>
  void AdaptorMaxOrder<ManagerImplementation>::update_self_helper() {
    static_assert(traits::MaxOrder > 2,
                  "No neighbourlist present; extension not possible.");

    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    this->nb_neigh.clear();
    this->offsets.clear();
    this->neighbours_atom_tag.clear();

    // #BUG8486@(markus) I now append the ghost atoms to the cluster index
    // container
    for (auto atom : this->manager) {
      //  Order 1, but variable Order is at 0, atoms, index 0
      using AddOrderLoop = AddOrderLoop<atom.order(), traits::NeighbourListType,
                                        IsCompactCluster,
                                        atom.order() == (traits::MaxOrder - 1)>;

      auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};

      auto indices = atom.get_cluster_indices();
      atom_cluster_indices.push_back(indices);

      AddOrderLoop::loop(atom, atom, 0, *this);
    }
    // correct the offsets for the new cluster order
    this->set_offsets();

    // add correct cluster_indices for the highest order
    auto & max_cluster_indices{
        std::get<traits::MaxOrder - 1>(this->cluster_indices_container)};
    max_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Returns the linear indices of the clusters (whose atom tags are stored
   * in counters). For example when counters is just the list of atoms, it
   * returns the index of each atom. If counters is a list of pairs of indices
   * (i.e. specifying pairs), for each pair of indices i,j it returns the
   * number entries in the list of pairs before i,j appears.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  size_t AdaptorMaxOrder<ManagerImplementation>::get_offset_impl(
      const std::array<size_t, Order> & counters) const {
    static_assert(Order < traits::MaxOrder,
                  "this implementation handles only up to the respective"
                  " MaxOrder");
    // Order accessor: 0 - atoms
    //                 1 - pairs
    //                 2 - triplets
    //                 etc.
    // Order is determined by the ClusterRef building iterator, not by the
    // Order of the built iterator

    // necessary construct for static branching
    using IncreaseHelper_t =
        internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

    if (Order < (traits::MaxOrder - 1)) {
      // If not accessible at this order, call lower Order offsets from lower
      // order manager or push through to lower levels, if adaptors are
      // stacked.
      return IncreaseHelper_t::get_offset(*this->manager, counters);
    } else {
      // Counters is an array to call parent offset multiplet. This can then
      // be used to access the actual offset for the Order which was built
      // here. It needs to be cast into a smaller one to access the order of
      // this Cluster(Order-1) from the manager below
      std::array<size_t, Order - 1> counters_below{};
      for (size_t c_index{0}; c_index < Order - 1; ++c_index) {
        counters_below[c_index] = counters[c_index];
      }
      // Linear index of the Cluster (Order-1)
      auto i{this->manager->get_offset_impl(counters_below)};
      // Number of cluster in its current iteration
      auto j{counters[Order - 1]};
      auto tuple_index{i + j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
