/**
 * file   adaptor_full_neighbour_list.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   16 Oct 2018
 *
 * @brief implements an adaptor for structure_managers, extending and existing
 * half/minimal neighbour list to a full one, where both permutations of pairs
 * exist.
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

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"


#ifndef ADAPTOR_FULL_LIST_H
#define ADAPTOR_FULL_LIST_H

namespace rascal {
  /*
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorFullList;

  /*
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
   * This adaptor ensures guarantiees, that each pair is contained only once
   * without a permutation.
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
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using PairRef_t = ClusterRef_t<2>;

    static_assert(traits::MaxOrder > 1,
                  "AdaptorFullList needs pairs.");
    static_assert(traits::MaxOrder < 3,
                  "AdaptorFullList does not work with Order > 2.");

    //! Default constructor
    AdaptorFullList() = delete;

    /**
     * Extend a half neighbour list to a full neighbour list, where both
     * permutations of a pair are present.
     */
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
    template<class ... Args>
    void update(Args&&... arguments);

    /**
     * returns the cutoff from the underlying manager which built the
     * neighbourlist
     */
    inline double get_cutoff() const {return this->manager.get_cutoff();}

    // //! returns the distance between atoms in a given pair
    // template <size_t Order, size_t Layer>
    // inline const double & get_distance(const ClusterRefKey<Order, Layer> &
    //                                    pair) const {
    //   return this->distance[pair];
    // }

    // template <size_t Order, size_t Layer>
    // inline double & get_distance(const ClusterRefKey<Order, Layer>& pair) {
    //   return this->distance[pair];
    // }

    inline size_t get_nb_clusters(int cluster_size) const {
      return this->atom_indices[cluster_size-1].size();
    }

    inline size_t get_size() const {
      return this->get_nb_clusters(1);
    }

    inline Vector_ref get_position(const int & index) {
      auto && original_index{this->atom_indices[0][index]};
      return this->manager.get_position(original_index);
    }

    //! get atom_index of index-th neighbour of this cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
				     & cluster,
				     int index) const {
      static_assert(Order <= 2,
                    "This implementation only handles upto pairs.");
      auto && offset = this->offsets[Order][cluster.get_cluster_index(Layer)];
      return this->atom_indices[Order][offset + index];
    }

    //! get atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent & /*parent*/,
				     size_t index) const {
      return this->atom_indices[0][index];
    }

    //! return atom type
    inline int & get_atom_type(const AtomRef_t & atom) {
      /**
       * careful, atom refers to our local index, for the manager, we need its
       * index:
       */
      auto && original_atom{this->atom_indices[0][atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    //! return atom type, const ref
    inline const int & get_atom_type(const AtomRef_t & atom) const {
      /**
       * careful, atom refers to our local index, for the manager, we need its
       * index:
       */
      auto && original_atom{this->atom_indices[0][atom.get_index()]};
      return this->manager.get_atom_type(original_atom);
    }

    //! Returns atom type given an atom index
    inline int & get_atom_type(const int & atom_id) {
      auto && original_id{this->atom_indices[0][atom_id]};
      return this->manager.get_atom_type(original_id);
    }

    //! Returns a constant atom type given an atom index
    inline const int & get_atom_type(const int & atom_id) const {
      auto && original_id{this->atom_indices[0][atom_id]};
      return this->manager.get_atom_type(original_id);
    }

    /**
     * Returns the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
				  & counters) const {
      return this->offsets[Order][counters.back()];
    }

    //! Returns the number of neighbours of a given cluster
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
				   & cluster) const {
      static_assert(Order <= traits::MaxOrder-1,
                    "this implementation only handles atoms and pairs");
      return this->nb_neigh[Order][cluster.back()];
    }

  protected:
    /**
     * main function during construction of a the full neighbourlist.
     * @param atom the atom to add to the list
     * @param Order select whether it is an i-atom (order=1), j-atom (order=2),
     * or ...
     */
    template <size_t Order>
    inline void add_atom(int atom_index) {
      static_assert(Order <= traits::MaxOrder,
                    "you can only add neighbours to the n-th degree defined by "
                    "MaxOrder of the underlying manager");

      // add new atom at this Order
      this->atom_indices[Order].push_back(atom_index);
      // count that this atom is a new neighbour
      this->nb_neigh[Order].back()++;
      this->offsets[Order].back()++;

      for (auto i{Order+1}; i < traits::MaxOrder; ++i) {
        // make sure that this atom starts with zero lower-Order neighbours
        this->nb_neigh[i].push_back(0);
        // update the offsets
        this->offsets[i].push_back(this->offsets[i].back() +
                                   this->nb_neigh[i-1].back());
      }
    }

    template <size_t Order>
    inline void add_atom(const typename ManagerImplementation::template
                         ClusterRef<Order> & cluster) {
      return this->template add_atom <Order-1>(cluster.back());
    }

    // Reference to the underlying manager
    ManagerImplementation & manager;

    /**
     * store atom indices per order,i.e.
     *   - atom_indices[0] lists all i-atoms
     *   - etc
     */
    // TODO this can be hardcoded to MaxOrder=2?
    std::array<std::vector<int>, traits::MaxOrder> atom_indices;
    /**
     * store the number of j-atoms for every i-atom (nb_neigh[1])
     */
    std::array<std::vector<size_t>, traits::MaxOrder> nb_neigh;
    /**
     * store the offsets from where the nb_neigh can be counted
     */
    std::array<std::vector<size_t>, traits::MaxOrder>  offsets;
  private:
  };

  //----------------------------------------------------------------------------//
  template <class ManagerImplementation>
  AdaptorFullList<ManagerImplementation>::
  AdaptorFullList(ManagerImplementation & manager):
    manager{manager},
    atom_indices{},
    nb_neigh{},
    offsets{}
  {}

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorFullList<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void AdaptorFullList<ManagerImplementation>::update() {

    //! Reset cluster_indices for adaptor to fill with push back.
    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    //! initialise the neighbourlist
    for (size_t i{0}; i < traits::MaxOrder; ++i) {
      this->atom_indices[i].clear();
      this->nb_neigh[i].resize(0);
      this->offsets[i].resize(0);
    }
    this->nb_neigh[0].push_back(0);
    for (auto & vector: this->offsets) {
      vector.push_back(0);
    }

    // fill the list, at least pairs are mandatory for this to work
    auto & atom_cluster_indices{std::get<0>(this->cluster_indices_container)};
    auto & pair_cluster_indices{std::get<1>(this->cluster_indices_container)};

    size_t pair_counter{0};
    for (auto atom: this->manager) {
      this->add_atom(atom);
      /**
       * Add new depth layer for atoms (see LayerByOrder for
       * possible optimisation).
       */

      constexpr auto AtomLayer{
        compute_cluster_layer<atom.order()>
          (typename traits::LayerByOrder{})
          };

      Eigen::Matrix<size_t, AtomLayer+1, 1> indices;
      indices.template head<AtomLayer>() = atom.get_cluster_indices();
      indices(AtomLayer) = indices(AtomLayer-1);
      atom_cluster_indices.push_back(indices);

      for (auto pair: atom) {
        constexpr auto PairLayer{
          compute_cluster_layer<pair.order()>
            (typename traits::LayerByOrder{})};
        // add first pair anyways
        this->add_atom(pair);

        Eigen::Matrix<size_t, PairLayer+1, 1> indices_pair;
        indices_pair.template head<PairLayer>() = pair.get_cluster_indices();
        indices_pair(PairLayer) = pair_counter;
        pair_cluster_indices.push_back(indices_pair);

        pair_counter++;

        //! TODO: access the second index's pairs and add this one to it. It is
        // probably better to first collect all pairs and then write them
        // contiguously

        //! Now create the second pair by permutation of the
        // existing indices auto index_i{atom.back()}; auto
        // index_j{pair.back()};
      }
    }
  }
}  // rascal

#endif /* ADAPTOR_FULL_LIST_H */
