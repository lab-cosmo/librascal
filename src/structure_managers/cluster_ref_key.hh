/**
 * @file   cluster_ref_key.hh
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

#ifndef SRC_STRUCTURE_MANAGERS_CLUSTER_REF_KEY_HH_
#define SRC_STRUCTURE_MANAGERS_CLUSTER_REF_KEY_HH_

#include "structure_managers/cluster_ref_base.hh"

#include <Eigen/Dense>

#include <array>
#include <iostream>
#include <tuple>

namespace rascal {
  /* ---------------------------------------------------------------------- */
  /**
   * Utilities for layer calculations and manipulations.
   */

  /**
   * Increases each layer of all existing orders.
   *
   * @tparam T the std::index_sequence representing the layer by order sequence.
   */
  template <class T>
  struct LayerIncreaser {};

  /**
   * @tparam LayersByOrder the layer by orders as variadic arguments
   */
  template <size_t... LayersByOrder>
  struct LayerIncreaser<std::index_sequence<LayersByOrder...>> {
    using type = std::index_sequence<(LayersByOrder + 1)...>;
  };

  template <size_t... LayersByOrder>
  using LayerIncreaser_t =
      typename LayerIncreaser<std::index_sequence<LayersByOrder...>>::type;

  /* ---------------------------------------------------------------------- */

  /**
   * Extends the layer by order index sequence by an additional cluster order.
   * Used when order is increased by an adaptor. See traits of`AdaptorMaxOrder`.
   *
   * @tparam T the std::index_sequence representing the layer by order sequence.
   */
  template <class T>
  struct LayerExtender;

  /**
   * @tparam LayersByOrder the layer by orders as variadic arguments
   */
  template <size_t... LayersByOrder>
  struct LayerExtender<std::index_sequence<LayersByOrder...>> {
    using type = std::index_sequence<LayersByOrder..., 0>;
  };

  template <size_t... LayersByOrder>
  using LayerExtender_t =
      typename LayerExtender<std::index_sequence<LayersByOrder...>>::type;

  /* ---------------------------------------------------------------------- */

  /**
   * Static access at position `Index` in the type of the parameter
   * `index_sequence`.
   *
   * @tparam Index
   * @param index_sequence
   *
   * @return returns the `Index` elment in the index sequence.
   */
  template <size_t Index, size_t... Ints>
  constexpr size_t get_index_from_sequence(
      const std::index_sequence<Ints...> & /*index_sequence*/) {
    return std::get<Index>(std::make_tuple(Ints...));
  }

  /**
   * From the index sequence given by the type of `layers_by_order` the layer is
   * returned at order/position `Order`.
   *
   * @tparam Order
   * @param layers_by_order the layers of each order in an index sequence
   *
   * @return returns the layer at the order `Order`
   */
  template <size_t Order, size_t... Ints>
  constexpr size_t
  get_layer(const std::index_sequence<Ints...> & layers_by_order) {
    static_assert(Order > 0, "Order is <1 this should not be.");
    static_assert(Order <= sizeof...(Ints),
                  "Order should be within the MaxOrder.");
    return get_index_from_sequence<Order - 1>(layers_by_order);
  }

  /**
   * Returns the index sequence given by the `Ints` as array.
   *
   * @param layers_by_order the type captures the layers of each order in an
   *        index sequence.
   *
   * @return returns `Ints` as array
   */
  template <size_t... Ints>
  constexpr std::array<size_t, sizeof...(Ints)> index_sequence_to_array(
      const std::index_sequence<Ints...> & /*layers_by_order*/) {
    return std::array<size_t, sizeof...(Ints)>{{Ints...}};
  }

  /**
   * Returns the smallest number of layers up to the order `ActiveMaxOrder` in
   * the `LayersByOrder` index sequence. It is needed to access properties at a
   * specific layer in a manager stack.
   *
   * @tparam ActiveMaxOrder the maximum order to check the `LayersByOrders` for
   *         the minimum value.
   * @param layers_by_order the type captures the layers of each order in an
   *        index sequence.
   *
   * @return the smallest value in the given index sequence up to Order elements
   */
  template <size_t ActiveMaxOrder, size_t... LayersByOrder>
  constexpr size_t get_min_layer(
      const std::index_sequence<LayersByOrder...> & /*layers_by_order*/) {
    static_assert(ActiveMaxOrder > 0,
                  "ActiveMaxOrder should be greater than 0.");
    static_assert(ActiveMaxOrder <= sizeof...(LayersByOrder),
                  "ActiveMaxOrder should not greater than the MaxOrder.");
    // transforms the LayersByOrder to an array
    //constexpr size_t arr[] = {LayersByOrder...};
    constexpr std::array<size_t, sizeof...(LayersByOrder)> arr {{
        LayersByOrder...}};
    return *std::min_element(std::begin(arr), std::begin(arr) + ActiveMaxOrder);
  }

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
   * For these reasons ClusterRefKey is templated by two arguments: Order that
   * specifies the number of atoms in the cluster, and Layer that specifies how
   * many layers of managers/adaptors are stacked at the point at which the
   * cluster reference is introduced.
   */
  template <size_t Order, size_t Layer>
  class ClusterRefKey : public ClusterRefBase {
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
    using AtomIndex_t = std::array<int, Order>;

    //! Default constructor
    ClusterRefKey() = delete;

    /**
     * direct constructor. Initialized with an array of atoms indices,
     * and a cluster reference data
     */
    ClusterRefKey(AtomIndex_t atom_tag_list, IndexConstArray cluster_indices)
        : Parent{Order, Layer}, atom_tag_list{atom_tag_list},
          cluster_indices{cluster_indices.data()} {}

    //! Copy constructor
    ClusterRefKey(const ClusterRefKey & other) = default;

    //! Move constructor
    ClusterRefKey(ClusterRefKey && other) = default;

    //! Destructor
    virtual ~ClusterRefKey() = default;

    //! Copy assignment operator
    ClusterRefKey & operator=(const ClusterRefKey & other) = delete;

    //! Move assignment operator
    ClusterRefKey & operator=(ClusterRefKey && other) = default;

    //! returns the atom tags of the current cluster
    const AtomIndex_t & get_atom_tag_list() const {
      return this->atom_tag_list;
    }

    //! returns the first atom tag in this cluster
    int front() const { return this->atom_tag_list.front(); }
    //! returns the last atom tag in this cluster
    int back() const { return this->atom_tag_list.back(); }
    /* the internal cluster neighbour is the neighbour which was added as
     * neighbour in the creation of this cluster
     */
    int get_internal_neighbour_atom_tag() const { return this->back(); }

    /*
     * From an cluster of form (i,j,..., n) it returns the tag of atom n
     */
    int get_atom_tag() const { return this->back(); }

    //! returns the cluster's index, given a specific layer
    size_t get_cluster_index(const size_t layer) const {
      return this->cluster_indices(layer);
    }

    //! returns the complete cluster indices (stacking history)
    IndexConstArray get_cluster_indices() const {
      return this->cluster_indices;
    }

    //! returns the order of the current cluster
    constexpr static size_t order() { return Order; }

    //! returns the layer of the current cluster
    constexpr static size_t cluster_layer() { return Layer; }

   protected:
    /**
     *  Array with unique atom tags. These can be user defined to refer to
     *  the exact same atom, e.g. in a Monte-Carlo simulation, where atoms are
     *  swapped.
     */
    AtomIndex_t atom_tag_list;
    /**
     * Cluster indices by layer order, highest layer, means last adaptor, and
     * means last entry (.back())
     */
    IndexConstArray cluster_indices;
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_CLUSTER_REF_KEY_HH_
