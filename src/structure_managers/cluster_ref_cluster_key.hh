/**
 * file   cluster_ref_cluster_key.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricer <markus.stricker@epfl.ch>
 *
 * @date   21 Jun 2018
 *
 * @brief an accessor class for getting access to clusters along a stack of
 *        neighbourhood/adaptors
 *
 * Copyright  2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_CLUSTER_REF_CLUSTER_KEY_HH_
#define SRC_STRUCTURE_MANAGERS_CLUSTER_REF_CLUSTER_KEY_HH_

#include "structure_managers/cluster_ref_base.hh"

#include <Eigen/Dense>

#include <tuple>
#include <array>
#include <iostream>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  /**
   * Accessor class for a reference to a cluster, i.e. a tuple of atoms (atoms,
   * pairs, triples, ...). The reference does not store data about the actual
   * tuple, just contains all the information needed to locate the infor in the
   * appropriate arrays that are stored in a Manager class.
   *
   * Given that Manager classes can be 'stacked', e.g. using a strict cutoff on
   * top of a loose neighbor list, the reference must know in which order of the
   * hierarchy the data.
   *
   * For these reasons ClusterRefClusterKey is templated by two arguments: Order that
   * specifies the number of atoms in the cluster, and Layer that specifies how
   * many layers of managers/adaptors are stacked at the point at which the
   * cluster reference is introduced.
   */
  template <size_t Order, size_t Layer>
  class ClusterRefClusterKey : public ClusterRefBase {
   public:
    /**
     * Index array types need both a constant and a non-constant version. The
     * non-const version can and needs to be cast into a const version in
     * argument.
     */
    using Parent = ClusterRefBase;
    using IndexConstArray =
        Eigen::Map<const Eigen::Matrix<size_t, Layer + 1, 1>>;
    using IndexArray = Eigen::Map<Eigen::Matrix<size_t, Layer + 1, 1>>;
    // TODO(alex) it is used as size_t to store atom indices in AMO, change
    // this to be consistent and change this here to size_t because indices can
    // not be negative
    using AtomIndex_t = std::array<int, Order>;

    //! Default constructor
    ClusterRefClusterKey() = delete;

    /**
     * direct constructor. Initialized with an array of atoms indices,
     * and a cluster reference data
     */
    ClusterRefClusterKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices)
        : Parent{Order, Layer}, atom_indices{atom_indices},
          cluster_indices{cluster_indices.data()} {}

    //! Copy constructor
    ClusterRefClusterKey(const ClusterRefClusterKey & other) = default;

    //! Move constructor
    ClusterRefClusterKey(ClusterRefClusterKey && other) = default;

    //! Destructor
    virtual ~ClusterRefClusterKey() = default;

    //! Copy assignment operator
    ClusterRefClusterKey & operator=(const ClusterRefClusterKey & other) = delete;

    //! Move assignment operator
    ClusterRefClusterKey & operator=(ClusterRefClusterKey && other) = default;

    //! returns the atom indices of the current cluster
    const inline AtomIndex_t & get_atom_indices() const {
      return this->atom_indices;
    }

    //! returns the first atom index in this cluster
    const int & front() const { return this->atom_indices.front(); }
    //! returns the last atom index in this cluster
    const int & back() const { return this->atom_indices.back(); }
    /* the internal cluster neighbour is the neighbour which was added as
     * neighbour in the creation of this cluster
     */
    // TODO(alex) rename get neighbour_atom_index
    const int & get_internal_cluster_neighbour_index() const {
      return this->back();
    }
    
    template<size_t Order_=Order>
    std::enable_if_t<Order_==1,const int &> get_atom_index() const {
      return this->back();
    }

    //! returns the cluster's index, given a specific layer
    inline size_t get_cluster_index(const size_t layer) const {
      return this->cluster_indices(layer);
    }

    //! returns the complete cluster indices (stacking history)
    inline IndexConstArray get_cluster_indices() const {
      return this->cluster_indices;
    }

    //! returns the order of the current cluster
    constexpr static inline size_t order() { return Order; }

    //! returns the layer of the current cluster
    constexpr static inline size_t cluster_layer() { return Layer; }

   protected:
    /**
     *  Array with unique atom indices. These can be user defined to refer to
     *  the exact same atom, e.g. in a Monte-Carlo simulation, where atoms are
     *  swapped.
     */
    AtomIndex_t atom_indices;
    /**
     * Cluster indices by layer order, highest layer, means last adaptor, and
     * means last entry (.back())
     */
    IndexConstArray cluster_indices;
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_CLUSTER_REF_CLUSTER_KEY_HH_
