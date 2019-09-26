/**
 * file   adaptor_full_neighbour_list.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   25 Oct 2018
 *
 * @brief implements an adaptor for structure_managers, extending the original
 *        manager so that the neighbourlist contains each pair twice, i.e. all
 *        permutations are present
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_FULL_NEIGHBOUR_LIST_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_FULL_NEIGHBOUR_LIST_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorFullList;

  /**
   * specialisation of traits for full neighbour list adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorFullList<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{
        ManagerImplementation::traits::Strict};
    constexpr static bool HasDistances{
        ManagerImplementation::traits::HasDistances};
    constexpr static bool HasDirectionVectors{
        ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};
    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
        AdaptorTraits::NeighbourListType::full};
    // New pairs are added at this layer, which did not exist before. Therefore
    // the layering has to be reset.
    constexpr static size_t AtomLayer{
        get<0>(typename LayerIncreaser<
               MaxOrder,
               typename ManagerImplementation::traits::LayerByOrder>::type{})};
    using LayerByOrder = std::index_sequence<AtomLayer, 0>;
  };

  /**
   * This adaptor guarantees, that each pair is contained twice, i.e. including
   * permutations.
   *
   * This interface should be implemented by all managers with the trait
   * AdaptorTraits::NeighbourListType{AdaptorTraits::NeighbourListType::full}
   */
  template <class ManagerImplementation>
  class AdaptorFullList
      : public StructureManager<AdaptorFullList<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorFullList<ManagerImplementation>> {
   public:
    using Parent = StructureManager<AdaptorFullList<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorFullList>;
    using Manager_t = AdaptorFullList<ManagerImplementation>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using parent_traits = typename ManagerImplementation::traits;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    // The stacking of this Adaptor is only possible on a manager which has a
    // pair list (MaxOrder=2). This is ensured here.
    static_assert(traits::MaxOrder > 1, "AdaptorFullList needs pairs.");
    static_assert(traits::MaxOrder < 3,
                  "AdaptorFullList does not work with Order > 2.");
    // TODO(markus): add this trait to all structure managers
    // static_assert(parent_traits::NeighbourListType
    //               == AdaptorTraits::NeighbourListType::half,
    //               "extends only a minimal neighbour list.");

    //! Default constructor
    AdaptorFullList() = delete;

    //! Extend a minimal/half neighbour list to a full neighbour list.
    explicit AdaptorFullList(ImplementationPtr_t manager);

    AdaptorFullList(ImplementationPtr_t manager, std::tuple<>)
        : AdaptorFullList(manager) {}

    AdaptorFullList(ImplementationPtr_t manager,
                    const Hypers_t & /*adaptor_hypers*/)
        : AdaptorFullList(manager) {}

    //! Copy constructor
    AdaptorFullList(const AdaptorFullList & other) = delete;

    //! Move constructor
    AdaptorFullList(AdaptorFullList && other) = default;

    //! Destructor
    virtual ~AdaptorFullList() = default;

    //! Copy assignment operator
    AdaptorFullList & operator=(const AdaptorFullList & other) = delete;

    //! Move assignment operator
    AdaptorFullList & operator=(AdaptorFullList && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update_self();

    //! update the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    /**
     * returns the cutoff from the underlying manager which built the
     * neighbourlist
     */
    inline double get_cutoff() const { return this->manager->get_cutoff(); }

    //! returns the number of atoms or pairs
    inline size_t get_nb_clusters(int order) const {
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

    //! returns the number of atoms
    inline size_t get_size() const { return this->get_nb_clusters(1); }

    //! returns the number of atoms
    inline size_t get_size_with_ghosts() const {
      return this->manager->get_size_with_ghosts();
    }

    //! returns position of the given atom tag
    inline Vector_ref get_position(const int & index) {
      return this->manager->get_position(index);
    }

    //! returns position of the given atom object
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager->get_position(atom.get_index());
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template <size_t Order, size_t Layer>
    inline int
    get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & cluster,
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
    inline int get_neighbour_atom_tag(const Parent &, size_t index) const {
      return this->manager->get_neighbour_atom_tag(*this->manager, index);
    }

    size_t get_atom_index(const int atom_tag) const {
      return this->manager->get_atom_index(atom_tag);
    }

    //! return atom type
    inline int & get_atom_type(const AtomRef_t & atom) {
      return this->manager->get_atom_type(atom.get_index());
    }

    //! return atom type, const ref
    inline const int & get_atom_type(const AtomRef_t & atom) const {
      return this->manager->get_atom_type(atom.get_index());
    }

    //! Returns atom type given an atom tag
    inline int & get_atom_type(const int & atom_id) {
      return this->manager->get_atom_type(atom_id);
    }

    //! Returns a constant atom type given an atom tag
    inline const int & get_atom_type(const int & atom_id) const {
      return this->manager->get_atom_type(atom_id);
    }

    /**
     * Returns the linear index of cluster (i.e., the count at which this
     * cluster appears in an iteration
     */
    template <size_t Order>
    inline size_t
    get_offset_impl(const std::array<size_t, Order> & counters) const {
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
    inline typename std::enable_if_t<(Order < (traits::MaxOrder - 1)), size_t>
    get_cluster_size_impl(const ClusterRefKey<Order, Layer> & cluster) const {
      return this->manager->get_cluster_size(cluster);
    }

    template <size_t Order, size_t Layer>
    inline typename std::enable_if_t<not(Order < traits::MaxOrder - 1), size_t>
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
      auto access_index = cluster.get_cluster_index(Layer);
      return nb_neigh[access_index];
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager() {
      return this->manager->get_shared_ptr();
    }

   protected:
    /* ---------------------------------------------------------------------- */
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
  AdaptorFullList<ManagerImplementation>::AdaptorFullList(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, nb_neigh{},
        neighbours_atom_tag{}, offsets{} {}

  /* ---------------------------------------------------------------------- */
  //! update function, which updates based on underlying manager
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorFullList<ManagerImplementation>::update(Args &&... arguments) {
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Update functions, which involes the extension of the neighbour list to one
   * which does includes all permutations of the pair
   */
  template <class ManagerImplementation>
  void AdaptorFullList<ManagerImplementation>::update_self() {
    // vector to locally gather all neighbours of an atom before building the
    // neighbour list
    std::vector<std::vector<int>> new_neighbours;

    // Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    // initialise empty data structures for the reduced neighbour list before
    // filling it
    this->nb_neigh.resize(0);
    this->offsets.resize(0);
    this->neighbours_atom_tag.resize(0);

    // prepare data structure to collect neighbours
    auto natoms = this->manager->get_size();
    new_neighbours.resize(natoms);
    for (auto & vector : new_neighbours) {
      // start with an empty list per atom
      vector.resize(0);
    }

    /* ---------------------------------------------------------------------- */
    // loop through all atoms and pairs and collect all neighbours in vector
    for (auto atom : *this->manager) {
      auto atom_tag{atom.get_atom_tag()};

      for (auto pair : atom) {
        auto neighbour_atom_index{
            this->get_atom_index(pair.get_internal_neighbour_atom_tag())};

        // add indices to their reciprocal list
        // -> already exists: this->new_neighbours[index_1].push_back(index_2);
        new_neighbours[neighbour_atom_index].push_back(atom_tag);
      }
    }

    /* ---------------------------------------------------------------------- */
    // reference to cluster indices
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    // build new neighbour list
    int offset{0};
    int pair_counter{0};

    for (auto atom : *this->manager) {
      auto atom_index = this->get_atom_index(atom.get_atom_tag());

      // Add new depth layer for atoms
      constexpr auto AtomLayer{
          compute_cluster_layer<atom.order()>(typename traits::LayerByOrder{})};

      Eigen::Matrix<size_t, AtomLayer + 1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer - 1);
      atom_cluster_indices.push_back(indices);

      int nneigh{0};
      for (auto pair : atom) {
        // add existing pairs
        auto neighbour_atom_tag = pair.get_internal_neighbour_atom_tag();
        this->neighbours_atom_tag.push_back(neighbour_atom_tag);
        nneigh++;

        // The layer of pairs is reinitialized with this adaptor. Therefore the
        // sice of the cluster indices is just "1". No need to copy underlying
        // indices, because they do not make sense in the stack.
        constexpr auto PairLayer{0};
        Eigen::Matrix<size_t, PairLayer + 1, 1> indices_pair;
        indices_pair(PairLayer) = pair_counter;
        pair_cluster_indices.push_back(indices_pair);
        pair_counter++;
      }

      // static expression for template parameter in cluster layer computation
      constexpr static auto PairOrder{2};
      // statically compute stacking height of pairs, which is to be increased
      // through extending the neighbour list
      constexpr static auto ActiveLayer{
          compute_cluster_layer<PairOrder>(typename traits::LayerByOrder{})};

      for (auto neighbour_atom_tag : new_neighbours[atom_index]) {
        this->neighbours_atom_tag.push_back(neighbour_atom_tag);
        nneigh++;

        Eigen::Matrix<size_t, ActiveLayer + 1, 1> indices_pair;
        // set cluster indices of the new pair to zero, since it does not exist
        // at the lower levels
        for (size_t i{0}; i < ActiveLayer; ++i) {
          indices_pair(i) = 0;
        }
        indices_pair(ActiveLayer) = pair_counter;
        pair_cluster_indices.push_back(indices_pair);
        pair_counter++;
      }
      // adjust offsets for correct access
      this->nb_neigh.push_back(nneigh);
      this->offsets.push_back(offset);
      offset += nneigh;
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_FULL_NEIGHBOUR_LIST_HH_
