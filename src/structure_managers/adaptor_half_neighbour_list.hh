/**
 * file   adaptor_half_neighbour_list.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   04 Oct 2018
 *
 * @brief implements an adaptor for structure_managers, filtering the
 * original manager so that the neighbourlist contains each pair
 * only once i.e. a half neighbour list
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

#ifndef ADAPTOR_HALF_LIST_H
#define ADAPTOR_HALF_LIST_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"

namespace rascal {
  /*
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorHalfList;

  /*
   * specialisation of traits for half neighbour list adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorHalfList<ManagerImplementation>> {

    constexpr static AdaptorTraits::Strict Strict{
      ManagerImplementation::traits::Strict};
    constexpr static bool HasDistances{
      ManagerImplementation::traits::HasDistances};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};

    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
      AdaptorTraits::NeighbourListType::half};
    // TODO: Future optimisation: do not increase depth for atoms
    // (they are all kept anyways, so no duplication necessary).
    using LayerByOrder = typename
      LayerIncreaser<MaxOrder,
                     typename
                     ManagerImplementation::traits::LayerByOrder>::type;
  };

  /**
   * This adaptor guarantees, that each pair is contained only once
   * without a permutation.
   *
   * This interface should be implemented by all managers with the trait
   * AdaptorTraits::NeighbourListType{AdaptorTraits::NeighbourListType::half}
   */
  template <class ManagerImplementation>
  class AdaptorHalfList: public
  StructureManager<AdaptorHalfList<ManagerImplementation>>
  {
  public:
    using Parent =
      StructureManager<AdaptorHalfList<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorHalfList>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using PairRef_t = ClusterRef_t<2>;

    /**
     * The stacking of this Adaptor is only possible on a manager which has a
     * pair list (MaxOrder=2). This is ensured here.
     */
    static_assert(traits::MaxOrder > 1,
                  "AdaptorHalfList needs pairs.");
    static_assert(traits::MaxOrder < 3,
                  "AdaptorHalfList does not work with Order > 2.");

    //! Default constructor
    AdaptorHalfList() = delete;

    /**
     * Reduce a full neighbour list to a half neighbour list (sometimes also
     * called minimal).
     */
    AdaptorHalfList(ManagerImplementation & manager);

    //! Copy constructor
    AdaptorHalfList(const AdaptorHalfList & other) = delete;

    //! Move constructor
    AdaptorHalfList(AdaptorHalfList && other) = default;

    //! Destructor
    virtual ~AdaptorHalfList() = default;

    //! Copy assignment operator
    AdaptorHalfList& operator=(const AdaptorHalfList & other) = delete;

    //! Move assignment operator
    AdaptorHalfList& operator=(AdaptorHalfList && other) = default;

    //! update just the adaptor assuming the underlying manager was updated
    void update();

    //! update the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    /**
     * returns the cutoff from the underlying manager which built the
     * neighbourlist
     */
    inline double get_cutoff() const {return this->manager.get_cutoff();}

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
        throw std::runtime_error("Can only handle single atoms and pairs");
      }
      }
    }

    inline size_t get_size() const {
      return this->get_nb_clusters(1);
    }

    //! Returns position of the given atom index
    inline Vector_ref get_position(const int & index) {
      return this->manager.get_position(index);
    }

    //! Returns position of the given atom object (useful for users)
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager.get_position(atom.get_index());
    }

    //! return position vector of the last atom in the cluster
    template<size_t Order, size_t Layer>
    inline Vector_ref get_neighbour_position(const ClusterRefKey<Order,
                                             Layer> & cluster) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "Order too large, not available.");

      return this->get_position(cluster.back());
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
				     & cluster,
				     size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");

      using IncreaseHelper_t =
        internal::IncreaseHelper<Order == (traits::MaxOrder-1)>;

      if (Order < (traits::MaxOrder-1)) {
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
      /**
       * careful, atom refers to our local index, for the manager, we need its
       * index:
       */
      return this->manager.get_atom_type(atom.get_index());
    }

    //! return atom type, const ref
    inline const int & get_atom_type(const AtomRef_t & atom) const {
      /**
       * careful, atom refers to our local index, for the manager, we need its
       * index:
       */
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
     * Returns the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
				  & counters) const {
      /**
       * The static assert with <= is necessary, because the template parameter
       * ``Order`` is one Order higher than the MaxOrder at the current
       * level. The return type of this function is used to build the next Order
       * iteration.
       */
      static_assert(Order <= traits::MaxOrder,
                    "this implementation handles only up to the respective"
                    " MaxOrder");
      /**
       * Order accessor: 0 - atoms
       *                 1 - pairs
       *                 2 - triplets
       *                 etc.
       * Order is determined by the ClusterRef building iterator, not by the Order
       * of the built iterator
       */
      return this->offsets[counters.front()];
    }

    //! Returns the number of neighbours of a given cluster
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
				   & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      /**
       * The static assert with <= is necessary, because the template parameter
       * ``Order`` is one Order higher than the MaxOrder at the current
       * level. The return type of this function is used to build the next Order
       * iteration.
       */

      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");

      if (Order < (traits::MaxOrder-1)) {
        return this->manager.get_cluster_size(cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
        return nb_neigh[access_index];
      }
    }

  protected:

    // Reference to the underlying manager
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

  //----------------------------------------------------------------------------//
  template <class ManagerImplementation>
  AdaptorHalfList<ManagerImplementation>::
  AdaptorHalfList(ManagerImplementation & manager):
    manager{manager},
    // atom_indices{},
    // nb_neigh{},
    // offsets{},
    nb_neigh{},
    neighbours{},
    offsets{}
  {}

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorHalfList<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorHalfList<ManagerImplementation>::update() {

    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    //! TEST
    /**
     * initialise empty data structures for the reduced neighbour list before
     * filling it
     */
    this->nb_neigh.resize(0);
    this->offsets.resize(0);
    this->neighbours.resize(0);

    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    int offset{0};
    size_t pair_counter{0};

    for (auto atom: this->manager) {
      /**
       * Add new depth layer for atoms (see LayerByOrder for possible
       * optimisation).
       */
      constexpr auto AtomLayer{
        compute_cluster_layer<atom.order()>
          (typename traits::LayerByOrder{})};

      Eigen::Matrix<size_t, AtomLayer+1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer-1);
      atom_cluster_indices.push_back(indices);

      auto index_i{atom.get_atom_index()};

      //! neighbours per atom counter to correct for offsets
      int nneigh{0};

      for (auto pair: atom) {
        constexpr auto PairLayer{
          compute_cluster_layer<pair.order()>
            (typename traits::LayerByOrder{})};

        auto index_j{pair.get_atom_index()};

        /**
         * This is the actual check for the half neighbour list: only pairs with
         * higher index_j than index_i are used. It should in principle ensure a
         * minimal neighbour list.
         */
        if (index_i < index_j) {

          this->neighbours.push_back(index_j);
          nneigh++;

          Eigen::Matrix<size_t, PairLayer+1, 1> indices_pair;
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
}  // rascal

#endif /* ADAPTOR_HALF_LIST_H */
