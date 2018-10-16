/**
 * file   adaptor_increase_maxorder.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jun 2018
 *
 * @brief implements an adaptor for structure_managers, which
 * creates a full and half neighbourlist if there is none and
 * triplets/quadruplets, etc. if existent.
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


#ifndef ADAPTOR_MAXORDER_H
#define ADAPTOR_MAXORDER_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include "lattice.hh"
#include "basic_types.hh"

#include <typeinfo>
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

    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
      ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    //! New MaxOrder upon construction!
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder+1};
    //! New Layer
    //! TODO: Is this the correct way to initialize the increased order?
    using LayerByOrder =
      typename LayerExtender<MaxOrder,
                             typename
                             ManagerImplementation::traits::LayerByOrder>::type;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Adaptor that increases the MaxOrder of an existing StructureManager. This
   * means, if the manager does not have a neighbourlist, there is nothing this
   * adaptor can do (hint: use adaptor_neighbour_list before and stack this on
   * top), if it exists, triplets, quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder: public
  StructureManager<AdaptorMaxOrder<ManagerImplementation>>
  {
  public:
    using Base = StructureManager<AdaptorMaxOrder<ManagerImplementation>>;

    using Parent =
      StructureManager<AdaptorMaxOrder<ManagerImplementation>>;
    using traits = StructureManager_traits<AdaptorMaxOrder>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template<size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    // TODO change these type to take them from the ManagerImplementation
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                                   Eigen::Dynamic>>;

    static_assert(traits::MaxOrder > 2,
                  "ManagerImplementation needs at least a pair list.");

    //! Default constructor
    AdaptorMaxOrder() = delete;

    /**
     * Given at least a pair list, this adaptor creates the next Order
     * list. I.e. from pairs to triplets, triplets to quadruplet, etc. Not
     * cutoff is needed, the cutoff is implicitly given by the neighbourlist,
     * which was built
     */
    AdaptorMaxOrder(ManagerImplementation & manager);

    //! Copy constructor
    AdaptorMaxOrder(const AdaptorMaxOrder & other) = delete;

    //! Move constructor
    AdaptorMaxOrder(AdaptorMaxOrder && other) = default;

    //! Destructor
    virtual ~AdaptorMaxOrder() = default;

    //! Copy assignment operator
    AdaptorMaxOrder& operator=(const AdaptorMaxOrder &other) = delete;

    //! Move assignment operator
    AdaptorMaxOrder& operator=(AdaptorMaxOrder &&other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes making triplets, quadruplets,
     * etc. depending on the MaxOrder, pair list has to be present.
     */
    void update();

    //! Updates the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    /**
     * Returns the linear indices of the clusters (whose atom indices are stored
     * in counters). For example when counters is just the list of atoms, it
     * returns the index of each atom. If counters is a list of pairs of indices
     * (i.e. specifying pairs), for each pair of indices i,j it returns the
     * number entries in the list of pairs before i,j appears.
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
				  & counters) const;

    //! Returns the number of clusters of size cluster_size
    inline size_t get_nb_clusters(size_t cluster_size) const {
      switch (cluster_size) {
      case traits::MaxOrder: {
        return this->neighbours.size();
        break;
      }
      default:
        return this->manager.get_nb_clusters(cluster_size);
        break;
      }
    }

    //! Returns number of clusters of the original manager
    inline size_t get_size() const {
      return this->manager.get_size();
    }

    //! Returns position of an atom with index atom_index
    inline Vector_ref get_position(const size_t & atom_index) {
      return this->manager.get_position(atom_index);
    }

    //! Returns position of the given atom object (useful for users)
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager.get_position(atom.get_index());
    }
    
    template<size_t Order, size_t Layer>
    inline Vector_ref get_neighbour_position(const ClusterRefKey<Order, Layer>
                                             & cluster) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");

      return this->get_position(cluster.back());
    }

    /**
     * Returns the id of the index-th (neighbour) atom of the cluster that is
     * the full structure/atoms object, i.e. simply the id of the index-th atom
     */
    inline int get_cluster_neighbour(const Parent& /*parent*/,
				     size_t index) const {
      return this->manager.get_cluster_neighbour(this->manager, index);
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
				     & cluster,
				     size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");
      if (Order < traits::MaxOrder-1) {
	      return this->manager.get_cluster_neighbour(cluster, index);
      } else {
	      auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
	      return this->neighbours[offset + index];
      }
    }

    //! Returns the number of neighbors of a given cluster
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");
      if (Order < traits::MaxOrder-1) {
	      return this->manager.get_cluster_size(cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
	      return nb_neigh[access_index];
      }
    }

    //! Returns the number of neighbors of an atom with the given index
    inline size_t get_cluster_size(const int & atom_index) const {
      return this->manager.get_cluster_size(atom_index);
    }

  protected:
    /**
     * Main function during extension of neighbourlist/triplet list.
     *
     * @param atom The atom to add to the list. Because the MaxOrder is
     * increased by one in this adaptor, the Order=MaxOrder
     */
    inline void add_atom(const int atom_index) {
      //! adds new atom at this Order
      this->atom_indices.push_back(atom_index);
      //! increases the number of neighbours
      this->nb_neigh.back()++;
      //! increases the offsets
      this->offsets.back()++;

      /**
       * extends the list containing the number of neighbours with a new 0 entry
       * for the added atom
       */
      this->nb_neigh.push_back(0);

      /**
       * extends the list containing the offsets and sets it with the number of
       * neighbours plus the offsets of the last atom
       */
      this->offsets.push_back(this->offsets.back() +
                              this->nb_neigh.back());
    }

    //! Extends the list containing the number of neighbours with a 0
    inline void add_entry_number_of_neighbours() {
      this->nb_neigh.push_back(0);
    }

    //! Adds a given atom index as new cluster neighbour
    inline void add_neighbour_of_cluster(const int atom_index) {
      //! adds `atom_index` to neighbours
      this->neighbours.push_back(atom_index);
      //! increases the number of neighbours
      this->nb_neigh.back()++;
    }

    //! Sets the correct offsets for accessing neighbours
    inline void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      this->offsets.reserve(n_tuples);
      this->offsets.resize(1);
      // this->offsets.push_back(0);
      for (size_t i{0}; i < n_tuples; ++i) {
        this->offsets.emplace_back(this->offsets[i] + this->nb_neigh[i]);
      }
    }

    /**
     * Interface of the add_atom function that adds the last atom in a given
     * cluster
     */
    template <size_t Order>
    inline void add_atom(const typename ManagerImplementation::template
                         ClusterRef<Order> & cluster) {
      static_assert(Order <= traits::MaxOrder,
                    "Order too high, not possible to add atom");
      return this->add_atom(cluster.back());
    }

    ManagerImplementation & manager;

    template<size_t Order, bool IsDummy> struct AddOrderLoop;

    //! Stores atom indices of current Order
    std::vector<size_t> atom_indices{}; //akin to ilist[]

    //! Stores the number of neighbours for every traits::MaxOrder-1-*plets
    std::vector<size_t> nb_neigh{};

    //! Stores all neighbours of traits::MaxOrder-1-*plets
    std::vector<size_t> neighbours{};

    //! Stores the offsets of traits::MaxOrder-1-*plets for accessing
    //! `neighbours`, from where nb_neigh can be counted
    std::vector<size_t> offsets{};

    size_t cluster_counter{0};
    
  private:
  };

  /* ---------------------------------------------------------------------- */
  //!Construction of the next level manager
  template <class ManagerImplementation>
  AdaptorMaxOrder<ManagerImplementation>::
  AdaptorMaxOrder(ManagerImplementation & manager):
    manager{manager},
    atom_indices{},
    nb_neigh{},
    offsets{}
  {
    if (traits::MaxOrder < 2) {
      throw std::runtime_error("No pair list in manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void AdaptorMaxOrder<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <size_t Order, bool IsDummy>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};
    using ClusterRef_t = typename ManagerImplementation::template
      ClusterRef<Order>;

    using NextOrderLoop = AddOrderLoop<Order+1,
				       (Order+1 == OldMaxOrder)>;

    static void loop(ClusterRef_t & cluster,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {

      //! do nothing, if MaxOrder is not reached, except call the next order
      for (auto next_cluster : cluster) {

        auto & next_cluster_indices
        {std::get<Order>(manager.cluster_indices_container)};
        next_cluster_indices.push_back(next_cluster.get_cluster_indices());

	NextOrderLoop::loop(next_cluster, manager);
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
   * and full lists. E.g. this has to be adjusted to include both, the i- and the
   * j-atoms of each pair as an i-atom in a triplet (center).
   */
  template <class ManagerImplementation>
  template <size_t Order>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop<Order, true> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};

    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;

    using traits = typename AdaptorMaxOrder<ManagerImplementation>::traits;

    static void loop(ClusterRef_t & cluster,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {

      /**
       * get all i_atoms to find neighbours to extend the cluster to the next
       * order
       */
      auto i_atoms = cluster.get_atom_indices();

      /**
       * vector of existing i_atoms in `cluster` to avoid doubling of atoms in
       * final list
       */
      std::vector<size_t> current_i_atoms;

      /**
       * a set of new neighbours for the cluster, which will be added to extend
       * the cluster
       */
      std::set<size_t> current_j_atoms;

      //! access to underlying manager for access to atom pairs
      auto & manager_tmp{cluster.get_manager()};


      for (auto atom_index : i_atoms) {
        current_i_atoms.push_back(atom_index);
        size_t access_index = manager.get_cluster_neighbour(manager,
                                                            atom_index);

        //! build a shifted iterator to constuct a ClusterRef<1>
        auto iterator_at_position{manager_tmp.get_iterator_at(access_index)};

        /**
         * ClusterRef<1> as dereference from iterator to get pairs of the
         * i_atoms
         */
        auto && j_cluster{*iterator_at_position};

        /**
         * collect all possible neighbours of the cluster: collection of all
         * neighbours of current_i_atoms
         */
        for (auto pair : j_cluster) {
          auto j_add = pair.back();
          if (j_add > i_atoms.back()) {
            current_j_atoms.insert(j_add);
          }
        }
      }

      //! delete existing cluster atoms from list to build additional neighbours
      std::vector<size_t> atoms_to_add{};
      std::set_difference(current_j_atoms.begin(), current_j_atoms.end(),
                          current_i_atoms.begin(), current_i_atoms.end(),
                          std::inserter(atoms_to_add, atoms_to_add.begin()));

      manager.add_entry_number_of_neighbours();
      if (atoms_to_add.size() > 0) {
        for (auto j: atoms_to_add) {
          manager.add_neighbour_of_cluster(j);
        }
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  /**
   * This is the loop, which runs recursively goes to the maximum Order and then
   * increases it by one (i.e. pairs->triplets, triplets->quadruplets, etc.
   */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::update() {
    /**
     * Increase an existing pair/triplet/quadruplet list to a higher Order,
     * loops to highest Order with AddOrderLoop and then adds all the pairs of
     * the given cluster.
     */

    static_assert(traits::MaxOrder > 2,
                  "No neighbourlist present; extension not possible.");

    for (auto atom : this->manager) {
      //! Order 1, variable Order is at 0, atoms, index 0
      using AddOrderLoop = AddOrderLoop<atom.order(),
                                        atom.order() == traits::MaxOrder-1>;
      auto & atom_cluster_indices{std::get<0>
          (this->cluster_indices_container)};
      atom_cluster_indices.push_back(atom.get_cluster_indices());
      AddOrderLoop::loop(atom, *this);
    }

    //! correct the offsets for the new cluster order
    this->set_offsets();
    //! add correct cluster_indices for the highest order
    auto & max_cluster_indices
    {std::get<traits::MaxOrder-1>(this->cluster_indices_container)};
    max_cluster_indices.fill_sequence();
  }


  /* ---------------------------------------------------------------------- */
  /* Returns the linear indices of the clusters (whose atom indices
   * are stored in counters). For example when counters is just the list
   * of atoms, it returns the index of each atom. If counters is a list of pairs
   * of indices (i.e. specifying pairs), for each pair of indices i,j it returns
   * the number entries in the list of pairs before i,j appears.
   */
  template<class ManagerImplementation>
  template<size_t Order>
  inline size_t AdaptorMaxOrder<ManagerImplementation>::
  get_offset_impl(const std::array<size_t, Order> & counters) const {

    static_assert(Order < traits::MaxOrder,
                  "this implementation handles only up to "
                  "the respective MaxOrder");
    /**
     * Order accessor: 0 - atoms
     *                 1 - pairs
     *                 2 - triplets
     *                 etc.
     * Order is determined by the ClusterRef building iterator, not by the Order
     * of the built iterator
     */
    if (Order == traits::MaxOrder-1) {
      /**
       * Counters as an array to call parent offset multiplet. This can then be
       * used to access the actual offset for the Order which was built here.
       */
      auto i{this->manager.get_offset_impl(counters)};
      auto j{counters[Order-1]};
      auto tuple_index{i+j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    } else {
      /**
       * If not accessible at this order, call lower Order offsets from lower
       * order manager or push through to lower levels, if adaptors are stacked.
       */
      return this->manager.get_offset_impl(counters);
    }
  }
}  // rascal

#endif /* ADAPTOR_MAXORDER_H */

/**
 * TODO: The construction of triplets is fine, but they occur multiple times. We
 * probably need to check for increasing atomic index to get rid of
 * duplicates. But this is in general a design decision, if we want full/half
 * neighbour list and full/half/whatever triplets and quadruplets
 */
