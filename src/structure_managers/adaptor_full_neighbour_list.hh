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
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef ADAPTOR_FULL_LIST_H
#define ADAPTOR_FULL_LIST_H

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
    // TODO: Future optimisation: do not increase depth for atoms
    // (they are all kept anyways, so no duplication necessary).
    using LayerByOrder = typename
      LayerIncreaser<MaxOrder,
                     typename
                     ManagerImplementation::traits::LayerByOrder>::type;
  };

  /**
   * This adaptor guarantees, that each pair is contained twice, i.e. including
   * permutations.
   *
   * This interface should be implemented by all managers with the trait
   * AdaptorTraits::NeighbourListType{AdaptorTraits::NeighbourListType::full}
   */
  template <class ManagerImplementation>
  class AdaptorFullList: public
  StructureManager<AdaptorFullList<ManagerImplementation>>
  {
  public:
    using Parent =
      StructureManager<AdaptorFullList<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorFullList>;
    using parent_traits = typename ManagerImplementation::traits;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;

    // The stacking of this Adaptor is only possible on a manager which has a
    // pair list (MaxOrder=2). This is ensured here.
    static_assert(traits::MaxOrder > 1,
                  "AdaptorFullList needs pairs.");
    static_assert(traits::MaxOrder < 3,
                  "AdaptorFullList does not work with Order > 2.");
    // TODO: add this trait to all structure managers
    // static_assert(parent_traits::NeighbourListType
    //               == AdaptorTraits::NeighbourListType::half,
    //               "extends only a minimal neighbour list.");

    //! Default constructor
    AdaptorFullList() = delete;

     //! Extend a minimal/half neighbour list to a full neighbour list.
    AdaptorFullList(ManagerImplementation & manager);

    //! Copy constructor
    AdaptorFullList(const AdaptorFullList & other) = delete;

    //! Move constructor
    AdaptorFullList(AdaptorFullList && other) = default;

    //! Destructor
    virtual ~AdaptorFullList() = default;

    //! Copy assignment operator
    AdaptorFullList& operator=(const AdaptorFullList & other) = delete;

    //! Move assignment operator
    AdaptorFullList& operator=(AdaptorFullList && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update();

    //! update the underlying manager as well as the adaptor
    template<class... Args>
    void update(Args&&... arguments);

    /**
     * returns the cutoff from the underlying manager which built the
     * neighbourlist
     */
    inline double get_cutoff() const {return this->manager.get_cutoff();}

    //! returns the number of atoms or pairs
    inline size_t get_nb_clusters(int cluster_size) const {
      switch (cluster_size) {
      case 1: {
        return this->manager.get_nb_clusters(cluster_size);
        break;
      }
      case 2: {
        return this->neighbours.size();
        break;
      }
      default: {
        throw std::runtime_error("Can only handle single atoms and pairs.");
      }
      }
    }

    //! returns the number of atoms
    inline size_t get_size() const {return this->get_nb_clusters(1);}

    //! returns position of the given atom index
    inline Vector_ref get_position(const int & index) {
      return this->manager.get_position(index);
    }

    //! returns position of the given atom object
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager.get_position(atom.get_index());
    }


    //! Returns the id of the index-th neighbour atom of a given cluster
    template<size_t Order, size_t Layer>
    inline int
    get_cluster_neighbour(const ClusterRefKey<Order, Layer> & cluster,
                          size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");

      // necessary helper construct for static branching
      using IncreaseHelper_t =
        internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

      if (Order < (traits::MaxOrder - 1)) {
        return IncreaseHelper_t::get_cluster_neighbour(this->manager, cluster,
                                                       index);
      } else {
        auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
        return this->neighbours[offset + index];
      }
    }

    //! get atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent & /*parent*/,
                                     size_t index) const {
      return this->manager.get_cluster_neighbour(this->manager, index);
    }

    //! return atom type
    inline int & get_atom_type(const AtomRef_t & atom) {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! return atom type, const ref
    inline const int & get_atom_type(const AtomRef_t & atom) const {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Returns atom type given an atom index
    inline int & get_atom_type(const int & atom_id) {
      return this->manager.get_atom_type(atom_id);
    }

    //! Returns a constant atom type given an atom index
    inline const int & get_atom_type(const int & atom_id) const {
      return this->manager.get_atom_type(atom_id);
    }

    /**
     * Returns the linear index of cluster (i.e., the count at which this
     * cluster appears in an iteration
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
                                  & counters) const {

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
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & cluster) const {
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

      if (Order < (traits::MaxOrder-1)) {
        return this->manager.get_cluster_size(cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
        return nb_neigh[access_index];
      }
    }

  protected:
    /* ---------------------------------------------------------------------- */
    //! Reference to the underlying manager
    ManagerImplementation & manager;

    //! Stores the number of neighbours for every atom after sorting
    std::vector<size_t> nb_neigh;

    //! Stores all neighbours, i.e. atom indices in a list
    std::vector<size_t> neighbours;

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
  AdaptorFullList<ManagerImplementation>::
  AdaptorFullList(ManagerImplementation & manager):
    manager{manager},
    nb_neigh{},
    neighbours{},
    offsets{}
  {}

  /* ---------------------------------------------------------------------- */
  //! update function, which updates based on underlying manager
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorFullList<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Update functions, which involes the extension of the neighbour list to one
   * which does includes all permutations of the pair
   */
  template <class ManagerImplementation>
  void AdaptorFullList<ManagerImplementation>::update() {

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
    this->neighbours.resize(0);

    // prepare data structure to collect neighbours
    auto natoms = manager.get_size();
    new_neighbours.resize(natoms);
    for (auto & vector : new_neighbours) {
        // start with an empty list per atom
        vector.resize(0);
    }

    /* ---------------------------------------------------------------------- */
    // loop through all atoms and pairs and collect all neighbours in vector
    for (auto atom: this->manager) {
      auto index_1{atom.get_atom_index()};

      for (auto pair: atom) {
        auto index_2{pair.get_atom_index()};

        // add indices to their reciprocal list
        // -> already exists: this->new_neighbours[index_1].push_back(index_2);
        new_neighbours[index_2].push_back(index_1);
      }
    }

    /* ---------------------------------------------------------------------- */
    // reference to cluster indices
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    // build new neighbour list
    int offset{0};
    int pair_counter{0};

    for (auto atom : this->manager) {
      auto index_i = atom.get_atom_index();

      // Add new depth layer for atoms
      constexpr auto AtomLayer{
        compute_cluster_layer<atom.order()>
          (typename traits::LayerByOrder{})};

      Eigen::Matrix<size_t, AtomLayer+1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer-1);
      atom_cluster_indices.push_back(indices);

      int nneigh{0};
      for (auto pair : atom) {
        // add existing pairs
        auto index_j = pair.get_atom_index();
        this->neighbours.push_back(index_j);
        nneigh++;

        constexpr auto PairLayer{
          compute_cluster_layer<pair.order()>
            (typename traits::LayerByOrder{})};

        Eigen::Matrix<size_t, PairLayer+1, 1> indices_pair;
        indices_pair.template head<PairLayer>() = pair.get_cluster_indices();
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

      for (auto index_j : new_neighbours[index_i]) {
        this->neighbours.push_back(index_j);
        nneigh++;

        Eigen::Matrix<size_t, ActiveLayer+1, 1> indices_pair;
        // set cluster indices of the new pair to zero, since it does not exist
        // at the lower levels
        // TODO: not sure, this is right
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
}  // rascal

#endif /* ADAPTOR_FULL_LIST_H */
