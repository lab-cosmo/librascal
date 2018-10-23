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

#ifndef SPECIES_MANAGER_H
#define SPECIES_MANAGER_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"

#include <typeinfo>
#include <set>
#include <vector>

namespace rascal {
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
  template <class ManagerImplementation>
  class SpeciesManager
  {
  public:
    using Manager_t = StructureManager<ManagerImplementation>;
    using traits = StructureManager_traits<SpeciesManager>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template<size_t Order>
    using ClusterRef_t =
      typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using Positions_ref = Eigen::Map<Eigen::Matrix<double, traits::Dim,
                                                   Eigen::Dynamic>>;

    static_assert(traits::MaxOrder > 1,
                  "ManagerImplementation needs an atom list.");

    //! Default constructor
    SpeciesManager() = delete;

    SpeciesManager(ManagerImplementation& manager, double cutoff);

    //! Copy constructor
    SpeciesManager(const SpeciesManager &other) = delete;

    //! Move constructor
    SpeciesManager(SpeciesManager &&other) = default;

    //! Destructor
    virtual ~SpeciesManager() = default;

    //! Copy assignment operator
    SpeciesManager& operator=(const SpeciesManager &other) = delete;

    //! Move assignment operator
    SpeciesManager& operator=(SpeciesManager &&other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes building either the neighbour list or to
     * make triplets, quadruplets, etc. depending on the MaxOrder
     */
    void update();

    //! Updates the underlying manager as well as the adaptor
    template<class ... Args>
    void update(Args&&... arguments);

    //! Returns cutoff radius of the neighbourhood manager
    inline double get_cutoff() const {return this->cutoff;}

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
    //?
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
      static_assert(Order < traits::MaxOrder,
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
	return this->manager.get_cluster_neighbour(cluster, index);
    }

    //! Returns atom type given an atom object AtomRef
    inline int & get_atom_type(const AtomRef_t& atom) {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Returns a constant atom type given an atom object AtomRef
    inline const int & get_atom_type(const AtomRef_t& atom) const {
      return this->manager.get_atom_type(atom.get_index());
    }

    //! Returns the number of neighbors of a given cluster
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");
      return this->manager.get_cluster_size(cluster);
    }

    //! Returns the number of neighbors of an atom with the given index
    inline size_t get_cluster_size(const int & atom_index) const {
      return this->manager.get_cluster_size(atom_index);
    }

  protected:

    ManagerImplementation & manager;

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
  //! constructor
  template <class ManagerImplementation>
  SpeciesManager<ManagerImplementation>::
  SpeciesManager(ManagerImplementation & manager, double cutoff):
    manager{manager},
    atom_indices{},
    nb_neigh{},
    offsets{}
  {
    if (traits::MaxOrder < 1) {
      throw std::runtime_error("No atoms in manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <class ... Args>
  void SpeciesManager<ManagerImplementation>::update(Args&&... arguments) {
    this->manager.update(std::forward<Args>(arguments)...);
    this->update();
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  void SpeciesManager<ManagerImplementation>::update() {
    /**
     * Standard case, increase an existing neighbour list or triplet list to a
     * higher Order
     */

    //TODO
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
  inline size_t SpeciesManager<ManagerImplementation>::
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
    if (Order == 1) {
      return this->offsets[counters.front()];
    } else if (Order == traits::MaxOrder-1) {
      /**
       * Counters as an array to call parent offset multiplet. This can then be
       * used to access the actual offset for the next Order here.
       */
      auto i{this->manager.get_offset_impl(counters)};
      auto j{counters[Order-1]};
      auto tuple_index{i+j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    } else {
      /**
       * If not accessible at this order, call lower Order offsets from lower
       * order manager(s).
       */
      return this->manager.get_offset_impl(counters);
    }
  }
}  // rascal

#endif /* SPECIES_MANAGER_H */

/**
 * Questions
 * - How to ensure to know the depth of 'explosion' is known?
 * - How to formulate the existing depth?
 * - Combinatorics: sorting before or after pairs/triplets?
 * - method .get_elements (), returns a tuple of the current elements in the cluster
 */
