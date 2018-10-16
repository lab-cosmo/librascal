/**
 * file   structure_manager_lammps.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef STRUCTURE_MANAGER_LAMMPS_H
#define STRUCTURE_MANAGER_LAMMPS_H

#include "structure_managers/structure_manager.hh"

#include <stdexcept>
#include <vector>

namespace rascal {
  //! forward declaration for traits
  class StructureManagerLammps;

  /**
   * traits specialisation for Lammps manager The traits are used for vector
   * allocation and further down the processing chain to determine what
   * functionality the given StructureManager already contains to avoid
   * recomputation.  See also the implementation of adaptors.
   */
  template <>
  struct StructureManager_traits<StructureManagerLammps> {
    constexpr static int Dim{3};
    constexpr static size_t MaxOrder{2};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    using LayerByOrder = std::index_sequence<0, 0>;
  };

  //----------------------------------------------------------------------------//
  //! Definition of the new StructureManagerLammps class.
  class StructureManagerLammps:
    //! It inherits publicly everything from the base class
    public StructureManager<StructureManagerLammps>
  {
  public:
    using traits = StructureManager_traits<StructureManagerLammps>;
    using Parent = StructureManager<StructureManagerLammps>;
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;
    // template <size_t Order>
    // using ClusterRef_t = typename Parent::template ClusterRef<Order>;

    //! Default constructor
    StructureManagerLammps() = default;

    //! Copy constructor
    StructureManagerLammps(const StructureManagerLammps & other) = delete;

    //! Move constructor
    StructureManagerLammps(StructureManagerLammps && other) = default;

    //! Destructor
    virtual ~StructureManagerLammps() = default;

    //! Copy assignment operator
    StructureManagerLammps
    & operator=(const StructureManagerLammps & other) = delete;

    //! Move assignment operator
    StructureManagerLammps
    & operator=(StructureManagerLammps && other) = default;

    /**
     * resetting is required every time the list changes. Here, this
     * is implemented without explicit dependency to lammps. The
     * signature could be simplified by including lammps as a
     * dependency, but it is unclear that the convenience would
     * outweigh the hassle of maintaining the dependency.
     *
     * @param inum Property `inum` in the lammps `NeighList` structure
     *
     * @param tot_num sum of the properties `nlocal` and `nghost` in the
     *                lammps `Atom` structure
     *
     * @param ilist Property `ilist` in the lammps `NeighList` structure
     *
     * @param numneigh Property `numneigh` in the lammps `NeighList` structure
     *
     * @param firstneigh Property `firstneigh` in the lammps `NeighList` structure
     *
     * @param x Property `x` in the lammps `Atom` structure
     *
     * @param f Property `f` in the lammps `Atom` structure
     *
     * @param type Property `type` in the lammps `Atom` structure
     *
     * @param eatom per-atom energy
     *
     * @param vatom per-atom virial
     */
    void update(const int & inum, const int & tot_num,
                int * ilist, int * numneigh, int ** firstneigh,
                double ** x, double ** f, int * type,
                double * eatom, double ** vatom);



    //! return position vector of an atom given the atom index
    inline Vector_ref get_position(const size_t & atom_index) {
      auto * xval{this->x[atom_index]};
      return Vector_ref(xval);
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

    //! get atom type
    inline const int & get_atom_type(const int & atom_index) const {
      return this->type[atom_index];
    }

    //! Returns atom type given an atom index
    inline int & get_atom_type(const int & atom_index) {
      return this->type[atom_index];
    }

    //! return number of I atoms in the list
    inline size_t get_size() const {
      return this->inum;
    }

    //! return the number of neighbours of a given atom
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer> & cluster)
      const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      return this->numneigh[cluster.back()];
    }

    //! return the index-th neighbour of the last atom
    //! in a cluster with cluster_size = 1 (atoms)
    //! which can be used to construct pairs
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
                                     & cluster,
                                     size_t index) const {
      static_assert(Order == traits::MaxOrder-1,
                    "this implementation only handles atoms and identify its "
                    "index-th neighbour.");
      auto && i_atom_id{cluster.back()};
      return this->firstneigh[std::move(i_atom_id)][index];
    }

    /**
     * return the atom_index of the index-th atom in manager parent here is
     * dummy and is used for consistency in other words, atom_index is the
     * global LAMMPS atom index.
     */
    inline int get_cluster_neighbour(const Parent & /*cluster*/,
                                     size_t index) const {
      return this->ilist[index];
    }

    /**
     * provided an atom, returns the cumulative numbers of pairs
     * up to the first pair in which the atom is the I atom
     * this only works for atom
     */
    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
                                  & counters) const;

    /**
     * return the number of clusters of size cluster_size.
     * Can only handle cluster_size 1 (atoms) and cluster_size 2 (pairs).
    */
    size_t get_nb_clusters(int cluster_size) const;

  protected:
    int inum{};
    int tot_num{}; //includes ghosts
    int * ilist{};
    int * numneigh{};
    int ** firstneigh{};
    double **x{}; //! pointer to pointer
    double **f{};
    int * type{};
    double * eatom{};
    double ** vatom{};
    int nb_pairs{}; //! number of clusters with cluster_size=2 (pairs)
    std::vector<int> offsets{};

  private:
  };


  /*
   * provided an atom, returns the cumulative numbers of pairs
   * up to the first pair in which the atom is the I atom
   * this only works for atom
   */
  template<size_t Order>
  inline size_t StructureManagerLammps::
  get_offset_impl(const std::array<size_t, Order> & counters) const {
    static_assert (Order == 1, "this manager can only give the offset "
                   "(= starting index) for a pair iterator, given the i atom "
                   "of the pair");
      return this->offsets[counters.front()];
  }
}  // rascal

#endif /* STRUCTURE_MANAGER_LAMMPS_H */
