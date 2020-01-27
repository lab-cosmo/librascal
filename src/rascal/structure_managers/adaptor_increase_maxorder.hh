/**
 * @file   rascal/structure_managers/adaptor_increase_maxorder.hh
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

#ifndef SRC_RASCAL_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
#define SRC_RASCAL_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_

#include "rascal/structure_managers/lattice.hh"
#include "rascal/structure_managers/property.hh"
#include "rascal/structure_managers/structure_manager.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/utils/utils.hh"

#include <algorithm>
#include <set>
#include <vector>

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
        typename LayerExtender<typename parent_traits::LayerByOrder>::type;
    using PreviousManager_t = ManagerImplementation;
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
    using ConstImplementationPtr_t =
        const std::shared_ptr<const ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorMaxOrder>;
    using PreviousManager_t = typename traits::PreviousManager_t;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;
    //! Order added by the adaptor
    static constexpr size_t AdditionalOrder{traits::MaxOrder};

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

    /**
     * Returns the linear indices of the clusters (whose atom tags are stored
     * in counters). For example when counters is just the list of atoms, it
     * returns the index of each atom. If counters is a list of pairs of indices
     * (i.e. specifying pairs), for each pair of indices i,j it returns the
     * number entries in the list of pairs before i,j appears.
     */
    template <size_t Order, bool C = (Order < 2), std::enable_if_t<C, int> = 0>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const {
      return this->manager->get_offset(counters);
    }

    template <size_t Order, bool C = (Order < 2),
              std::enable_if_t<not(C), int> = 0>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const {
      auto j{counters.back()};
      auto main_offset{this->offsets[j]};
      return main_offset;
    }

    //! Returns the number of clusters of size cluster_size
    size_t get_nb_clusters(size_t order) const {
      /**
       * Note: The case for order=1 is abmiguous: one possible answer is the
       * number of centers the other possibility is the number of centers +
       * ghost atoms. Please use the get_size or get_size_with_ghosts member
       * functions
       */
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

    /**
     * Returns the id of the index-th neighbour atom of a given cluster.
     * This function is a helper to build ClusterRef of order 1 and 2.
     */
    template <size_t Layer>
    int get_neighbour_atom_tag(const ClusterRefKey<1, Layer> & cluster,
                               size_t index) const {
      return this->manager->get_neighbour_atom_tag(cluster, index);
    }

    /**
     * Return the atoms tags associated with the neighbors (triplet -> 2,
     * quadruplet -> 3) of cluster that have been added by this adaptor.
     * This function is a helper to build ClusterRef of order
     * AdditionalOrder.
     */
    template <size_t Layer>
    std::array<int, AdditionalOrder - 1>
    get_neighbour_atom_tag_current(const ClusterRefKey<1, Layer> & cluster,
                                   size_t index) const {
      auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
      return this->neighbours_atom_tag[offset + index];
    }

    size_t get_atom_index(const int atom_tag) const {
      return this->manager->get_atom_index(atom_tag);
    }

    /**
     * Returns the number of neighbors of a given atom at a given TargetOrder
     * use implementation of the previous manager when
     * TargetOrder != AddedOrder
     */
    template <size_t TargetOrder, size_t Order, size_t Layer>
        typename std::enable_if_t <
        TargetOrder<traits::MaxOrder, size_t> get_cluster_size_impl(
            const ClusterRefKey<Order, Layer> & cluster) const {
      return this->manager->template get_cluster_size<TargetOrder>(cluster);
    }

    //! Returns the number of neighbours of a given atom at a given TargetOrder
    //! at TargetOrder == MaxOrder use local data to give the nb of cluster
    template <size_t TargetOrder, size_t Order, size_t Layer>
    typename std::enable_if_t<TargetOrder == traits::MaxOrder, size_t>
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      constexpr auto current_layer{
          get_layer<TargetOrder>(typename traits::LayerByOrder{})};
      auto access_index = cluster.get_cluster_index(current_layer);
      return this->nb_neigh[access_index];
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager_impl() {
      return this->manager->get_shared_ptr();
    }

    //! Get the manager used to build the instance
    ConstImplementationPtr_t get_previous_manager_impl() const {
      return this->manager->get_shared_ptr();
    }

   protected:
    template <bool IsCompactCluster>
    void update_self_helper();

    //! Extends the list containing the number of neighbours with a 0
    void add_entry_number_of_neighbours() { this->nb_neigh.push_back(0); }

    //! Adds a given atom tag as new cluster neighbour
    void add_neighbour_of_cluster(
        const std::array<int, AdditionalOrder - 1> atom_tag) {
      // adds `atom_tag` to neighbours
      this->neighbours_atom_tag.push_back(atom_tag);
      // increases the number of neighbours
      this->nb_neigh.back()++;
    }

    //! Sets the correct offsets for accessing neighbours
    void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      if (n_tuples > 0) {
        this->offsets.emplace_back(0);
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
    struct ForwardClusterIndices;

    //! Stores the number of neighbours for every traits::MaxOrder-1-clusters
    std::vector<size_t> nb_neigh{};

    //! Stores all neighbours atom tag of traits::MaxOrder-1-clusters
    std::vector<std::array<int, AdditionalOrder - 1>> neighbours_atom_tag{};

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
  struct AdaptorMaxOrder<ManagerImplementation>::ForwardClusterIndices {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;

    using NextOrderLoop =
        ForwardClusterIndices<Order + 1, NeighbourListType, IsCompactCluster,
                              (Order + 1 <= OldMaxOrder)>;

    using AtomClusterRef_t =
        typename ManagerImplementation::template ClusterRef<1>;

    // do nothing, if MaxOrder is not reached, except call the next order
    static void loop(AtomClusterRef_t & atom,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {
      // copy the cluster indices of all orders below
      for (auto && cluster : atom.template get_clusters_of_order<Order>(0)) {
        auto & cluster_indices{
            std::get<Order - 1>(manager.cluster_indices_container)};
        // keep copying underlying cluster indices, they are not changed
        auto indices{cluster.get_cluster_indices()};
        cluster_indices.push_back(indices);
        NextOrderLoop::loop(atom, manager);
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * At desired MaxOrder (plus one), here is where the magic happens and the
   * neighbours of the same order are added as the Order+1.  add check for non
   * half neighbour list.
   */
  template <class ManagerImplementation>
  template <size_t Order, AdaptorTraits::NeighbourListType NeighbourListType,
            bool IsCompactCluster>
  struct AdaptorMaxOrder<ManagerImplementation>::ForwardClusterIndices<
      Order, NeighbourListType, IsCompactCluster, false> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};

    using AtomClusterRef_t =
        typename ManagerImplementation::template ClusterRef<1>;

    using traits = typename AdaptorMaxOrder<ManagerImplementation>::traits;

    //! loop through the orders to get to the maximum order, this is agnostic to
    //! the underlying MaxOrder, just goes to the maximum
    static void loop(AtomClusterRef_t & /*atom */,
                     AdaptorMaxOrder<ManagerImplementation> & /*manager */) {}
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
    using ForwardClusterIndices =
        ForwardClusterIndices<2, traits::NeighbourListType, IsCompactCluster,
                              true>;
    static constexpr bool HasCenterPairAndIsOrderTwo{traits::HasCenterPair and
                                                     traits::MaxOrder - 1 == 2};
    // size_t turns out to give 0 if false and 1 if true. When adding Order==3
    // to the manager the input manager could have center pairs and
    // ClusterStart avoid iterating over those
    static constexpr size_t ClusterStart{
        static_cast<size_t>(HasCenterPairAndIsOrderTwo)};

    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};

    this->nb_neigh.clear();
    this->offsets.clear();
    this->neighbours_atom_tag.clear();

    for (auto atom : this->manager) {
      // copy the underlying cluster indices of this->manager (not changed)
      // do order == 1
      auto indices{atom.get_cluster_indices()};
      atom_cluster_indices.push_back(indices);
      // do the rest
      ForwardClusterIndices::loop(atom, *this);

      this->add_entry_number_of_neighbours();
      std::array<int, traits::MaxOrder> new_tag_list{};
      // loop over the highest order available in this->manager
      for (auto cluster :
           atom.template get_clusters_of_order<traits::MaxOrder - 1>(
               ClusterStart)) {
        auto && tag_list{cluster.get_atom_tag_list()};
        // copy the tags from the previous order
        for (size_t ii{0}; ii < traits::MaxOrder - 1; ++ii) {
          new_tag_list[ii] = tag_list[ii];
        }
        for (auto pair : atom.template get_clusters_of_order<2>()) {
          // copy the atom tag of the new order to the last element of the new
          // tag list
          new_tag_list.back() = pair.get_atom_tag();
          if (traits::NeighbourListType ==
              AdaptorTraits::NeighbourListType::half) {
            if (IsCompactCluster) {
              throw std::runtime_error("Not implemented yet.");
            } else {
              // only add strictly lexicographicaly ordered tags, e.g.
              // (1,2,2) or (1,5,3) are not valid and (1,3,5) is valid
              if (new_tag_list[traits::MaxOrder - 1] >
                  new_tag_list[traits::MaxOrder - 2]) {
                std::array<int, traits::MaxOrder - 1> tags{};
                for (size_t ii{1}; ii < traits::MaxOrder; ++ii) {
                  tags[ii - 1] = new_tag_list[ii];
                }
                this->add_neighbour_of_cluster(tags);
              }
            }
          } else if (traits::NeighbourListType ==
                     AdaptorTraits::NeighbourListType::full) {
            if (IsCompactCluster) {
              throw std::runtime_error("Not implemented yet.");
            } else {
              // only add tags_list where none of the tags are equal to one
              // another, e.g. (3,1,2) is valid but (3,2,2) is not.
              std::array<int, traits::MaxOrder> new_tag_list_s{new_tag_list};
              std::sort(new_tag_list_s.begin(), new_tag_list_s.end());
              auto any_equal = std::adjacent_find(new_tag_list_s.begin(),
                                                  new_tag_list_s.end());
              if (any_equal == new_tag_list_s.end()) {
                std::array<int, traits::MaxOrder - 1> tags{};
                for (size_t ii{1}; ii < traits::MaxOrder; ++ii) {
                  tags[ii - 1] = new_tag_list[ii];
                }
                this->add_neighbour_of_cluster(tags);
              }
            }
          } else {
            throw std::runtime_error("Not implemented.");
          }
        }
      }
    }
    // correct the offsets for the new cluster order
    this->set_offsets();

    // add correct cluster_indices for the highest order
    auto & max_cluster_indices{
        std::get<traits::MaxOrder - 1>(this->cluster_indices_container)};
    max_cluster_indices.fill_sequence();
  }

}  // namespace rascal

#endif  // SRC_RASCAL_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
