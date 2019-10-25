/**
 * @file   adaptor_half_neighbour_list.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief implements an adaptor for structure_managers, filtering the original
 *        manager so that the neighbourlist contains each pair only once i.e. a
 *        half neighbour list
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_HALF_NEIGHBOUR_LIST_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_HALF_NEIGHBOUR_LIST_HH_

#include "rascal_utility.hh"
#include "structure_managers/property.hh"
#include "structure_managers/structure_manager.hh"

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorHalfList;

  /**
   * specialisation of traits for half neighbour list adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorHalfList<ManagerImplementation>> {
    using parent_traits = StructureManager_traits<ManagerImplementation>;
    constexpr static AdaptorTraits::Strict Strict{parent_traits::Strict};
    constexpr static bool HasDistances{parent_traits::HasDistances};
    constexpr static bool HasDirectionVectors{
        parent_traits::HasDirectionVectors};
    constexpr static bool HasCenterPair{parent_traits::HasCenterPair};
    constexpr static int StackLevel{parent_traits::StackLevel + 1};
    constexpr static int Dim{parent_traits::Dim};
    constexpr static size_t MaxOrder{parent_traits::MaxOrder};
    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
        AdaptorTraits::NeighbourListType::half};
    using LayerByOrder =
        typename LayerIncreaser<MaxOrder,
                                typename parent_traits::LayerByOrder>::type;
  };

  /**
   * This adaptor guarantees, that each pair is contained only once without
   * permutations.
   *
   * This interface should be implemented by all managers with the trait
   * AdaptorTraits::NeighbourListType{AdaptorTraits::NeighbourListType::half}
   */
  template <class ManagerImplementation>
  class AdaptorHalfList
      : public StructureManager<AdaptorHalfList<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorHalfList<ManagerImplementation>> {
   public:
    using Manager_t = AdaptorHalfList<ManagerImplementation>;
    using Parent = StructureManager<Manager_t>;
    using ManagerImplementation_t = ManagerImplementation;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorHalfList>;
    using parent_traits = typename ManagerImplementation::traits;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    // The stacking of this Adaptor is only possible on a manager which has a
    // pair list (MaxOrder=2). This is ensured here.
    static_assert(traits::MaxOrder > 1, "AdaptorHalfList needs pairs.");
    static_assert(traits::MaxOrder < 3,
                  "AdaptorHalfList does not work with Order > 2.");

    // TODO(markus): add this to all structure managers
    // static_assert(parent_traits::NeighbourListType
    //               == AdaptorTraits::NeighbourListType::full,
    //               "adapts only a full neighbour list.");

    //! Default constructor
    AdaptorHalfList() = delete;

    /**
     * Reduce a full neighbour list to a half neighbour list (sometimes also
     * called minimal).
     */
    explicit AdaptorHalfList(ImplementationPtr_t manager);

    AdaptorHalfList(ImplementationPtr_t manager, std::tuple<>)
        : AdaptorHalfList(manager) {}

    AdaptorHalfList(ImplementationPtr_t manager,
                    const Hypers_t & /*adaptor_hypers*/)
        : AdaptorHalfList(manager) {}

    //! Copy constructor
    AdaptorHalfList(const AdaptorHalfList & other) = delete;

    //! Move constructor
    AdaptorHalfList(AdaptorHalfList && other) = default;

    //! Destructor
    virtual ~AdaptorHalfList() = default;

    //! Copy assignment operator
    AdaptorHalfList & operator=(const AdaptorHalfList & other) = delete;

    //! Move assignment operator
    AdaptorHalfList & operator=(AdaptorHalfList && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update_self();

    //! update the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    /**
     * returns the cutoff from the underlying manager which built the
     * neighbourlist
     */
    double get_cutoff() const { return this->manager->get_cutoff(); }

    //! returns the number of atoms or pairs
    size_t get_nb_clusters(int order) const {
      switch (order) {
      case 1: {
        return this->manager->get_nb_clusters(order);
        break;
      }
      case 2: {
        return this->neighbours_atom_tag.size();
        break;
      }
      default: {
        throw std::runtime_error("Can only handle single atoms and pairs.");
      }
      }
    }

    //! returns the number of center atoms
    size_t get_size() const { return this->manager->get_size(); }

    //! returns the number of atoms
    size_t get_size_with_ghosts() const {
      return this->manager->get_size_with_ghosts();
    }

    //! returns position of the given atom tag
    Vector_ref get_position(int index) {
      return this->manager->get_position(index);
    }

    //! returns position of the given atom object
    Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager->get_position(atom.get_index());
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

    //! get atom_tag of the index-th atom in manager
    int get_neighbour_atom_tag(const Parent &, size_t index) const {
      return this->manager->get_neighbour_atom_tag(*this->manager, index);
    }

    size_t get_atom_index(const int atom_tag) const {
      return this->manager->get_atom_index(atom_tag);
    }

    //! return atom type, const ref
    int get_atom_type(const AtomRef_t & atom) const {
      return this->manager->get_atom_type(atom.get_index());
    }

    //! Returns a constant atom type given an atom tag
    int get_atom_type(int atom_id) const {
      return this->manager->get_atom_type(atom_id);
    }

    //! check whether neighbours of ghosts were considered
    bool get_consider_ghost_neighbours() const {
      return this->manager->get_consider_ghost_neighbours();
    }

    /**
     * Returns the linear index of cluster (i.e., the count at which this
     * cluster appears in an iteration
     */
    template <size_t Order>
    size_t get_offset_impl(const std::array<size_t, Order> & counters) const {
      // The static assert with <= is necessary, because the template parameter
      // ``Order`` is one Order higher than the MaxOrder at the current
      // level. The return type of this function is used to build the next Order
      // iteration.
      static_assert(Order <= traits::MaxOrder,
                    "this implementation handles only up to the respective"
                    " MaxOrder");

      // Order is determined by the ClusterRef building iterator, not by the
      // Order of the built iterator.
      return this->offsets[counters.front()];
    }

    //! Returns the number of neighbours of a given cluster
    template <size_t Order, size_t Layer>
    size_t
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      /*
       * The static assert with <= is necessary, because the template parameter
       * ``Order`` is one Order higher than the MaxOrder at the current
       * level. The return type of this function is used to build the next Order
       * iteration.
       */
      static_assert(Order <= traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");

      if (Order < (traits::MaxOrder - 1)) {
        return this->manager->get_cluster_size_impl(cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
        return nb_neigh[access_index];
      }
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager() {
      return this->manager->get_shared_ptr();
    }

   protected:
    //! Reference to the underlying manager
    ImplementationPtr_t manager;

    //! Stores the number of neighbours for every atom after sorting
    std::vector<size_t> nb_neigh;

    //! Stores all neighbours, i.e. atom tags in a list
    std::vector<int> neighbours_atom_tag;

    /**
     * Stores the offsets for accessing `neighbours`; this is the entry point in
     * ``neighbours`` for each atom, from where the number of neighbours
     * ``nb_neigh`` can be accessed
     */
    std::vector<size_t> offsets;

   private:
  };

  /* ---------------------------------------------------------------------- */
  //! constructor implementations
  template <class ManagerImplementation>
  AdaptorHalfList<ManagerImplementation>::AdaptorHalfList(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, nb_neigh{},
        neighbours_atom_tag{}, offsets{} {
    // this->manager->add_child(this->get_weak_ptr());
  }

  /* ---------------------------------------------------------------------- */
  //! update function, which updates based on underlying manager
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorHalfList<ManagerImplementation>::update(Args &&... arguments) {
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Update functions, which involes the reduction of the neighbour list to one
   * which does not include permutations of pairs
   */
  template <class ManagerImplementation>
  void AdaptorHalfList<ManagerImplementation>::update_self() {
    // Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    // initialise empty data structures for the reduced neighbour list before
    // filling it
    this->nb_neigh.resize(0);
    this->offsets.resize(0);
    this->neighbours_atom_tag.resize(0);

    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    int offset{0};
    // counter for total number of pairs (minimal list) for cluster_indices
    size_t pair_counter{0};

    for (auto atom : this->manager->with_ghosts()) {
      // Add new depth layer for atoms (see LayerByOrder for possible
      // optimisation).
      constexpr auto AtomLayer{
          compute_cluster_layer<atom.order()>(typename traits::LayerByOrder{})};

      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);

      auto index_i{atom.get_atom_tag()};

      // neighbours per atom counter to correct for offsets
      int nneigh{0};

      for (auto pair : atom) {
        constexpr auto PairLayer{compute_cluster_layer<pair.order()>(
            typename traits::LayerByOrder{})};

        auto index_j{pair.get_atom_tag()};

        // This is the actual check for the half neighbour list: only pairs with
        // higher index_j than index_i are used. It should in principle ensure a
        // minimal neighbour list.
        if (index_i < index_j) {
          this->neighbours_atom_tag.push_back(index_j);
          nneigh++;

          // get underlying cluster indices to add at this layer
          Eigen::Matrix<size_t, PairLayer + 1, 1> indices_pair;
          indices_pair.template head<PairLayer>() = pair.get_cluster_indices();
          indices_pair(PairLayer) = pair_counter;
          pair_cluster_indices.push_back(indices_pair);

          pair_counter++;
        }
      }
      this->nb_neigh.push_back(nneigh);
      this->offsets.push_back(offset);

      offset += nneigh;
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_HALF_NEIGHBOUR_LIST_HH_
