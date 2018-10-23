/**
 * file   cluster_ref_key.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   21 Jun 2018
 *
 * @brief an accessor class for getting access to clusters along a stack of
 *        neighbourhood/adaptors
 *
 * Copyright © 2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef CLUSTERREFKEY_H
#define CLUSTERREFKEY_H

#include "structure_managers/cluster_ref_base.hh"

#include <tuple>
#include <array>

#include <Eigen/Dense>

#include <iostream>

namespace rascal {
  /**
   * Layer calculations and manipulations
   */
  // Computes layer by cluster dimension for new adaptor layer
  template <size_t MaxOrder, class T>
  struct LayerIncreaser{};

  template <size_t MaxOrder, size_t... Ints>
  struct LayerIncreaser<MaxOrder,
                        std::index_sequence<Ints...>> {
    using type = std::index_sequence<(Ints+1)...>;
  };

  template <size_t MaxOrder, size_t... Ints>
  using LayerIncreaser_t =
    typename LayerIncreaser<MaxOrder, std::index_sequence<Ints...>>::type;

  /**
   * Extends layer by cluster for an additional cluster dimension
   */
  template <size_t MaxOrder, class T>
  struct LayerExtender{};

  template <size_t MaxOrder, size_t... Ints>
  struct LayerExtender<MaxOrder,
                       std::index_sequence<Ints...>>{
    using type = std::index_sequence<Ints..., 0>;
  };

  template <size_t MaxOrder, size_t... Ints>
  using LayerExtender_t =
    typename LayerExtender<MaxOrder, std::index_sequence<Ints...>>::type;

  /**
   * Dynamic access to all layers by cluster dimension (probably not
   * necessary)
   */
  template <size_t MaxOrder, size_t... Ints>
  constexpr std::array<size_t, MaxOrder>
  get_layers(std::index_sequence<Ints...>) {
    return std::array<size_t, MaxOrder>{Ints...};
  }

  namespace internal {
    template <size_t head, size_t... tail>
    struct Min {
      constexpr static
      size_t value{ head < Min<tail...>::value ? head : Min<tail...>::value};
    };

    template <size_t head>
    struct Min<head> {
      constexpr static size_t value{head};
    };

    template <class Sequence>
    struct MinExtractor {};

    template <size_t... Ints>
    struct MinExtractor<std::index_sequence<Ints...>> {
      constexpr static size_t value {Min<Ints...>::value};
    };

    template <size_t Order, class Sequence, size_t... Ints>
    struct HeadExtractor {};

    template <size_t... seq>
    struct HeadExtractorTail {
      using type = std::index_sequence<seq...>;
    };

    template <size_t Order, size_t head, size_t... tail, size_t... seq>
    struct HeadExtractor<Order, std::index_sequence<seq...>, head, tail...> {
      using Extractor_t = std::conditional_t
        <(Order > 1),
        HeadExtractor<Order-1,
                      std::index_sequence<seq..., head>,
                      tail...>,
        HeadExtractorTail<seq..., head>>;
      using type = typename Extractor_t::type;
    };
  }  // internal

  // TODO please explain the logic allowing to resolve the layer of a cluster
  template <size_t Order, size_t... Ints>
  constexpr size_t compute_cluster_layer(const std::index_sequence<Ints...> &) {
    using ActiveDimensions = typename internal::HeadExtractor
      <Order, std::index_sequence<>, Ints...>::type;
    return internal::MinExtractor<ActiveDimensions>::value;
  }

  /**
   * Dynamic access to layer by cluster dimension (possibly not
   * necessary)
   */
  template <size_t MaxOrder, size_t... Ints>
  constexpr size_t
  get_layer(size_t index, std::index_sequence<Ints...>) {
    constexpr size_t arr[] {Ints...};
    return arr [index];
  }

  /**
   * Static access to layer by cluster dimension (e.g., for defining
   * template parameter `NbRow` for a property
   */
  template <size_t index, size_t... Ints>
  constexpr size_t get(std::index_sequence<Ints...>) {
    return get<index>(std::make_tuple(Ints...));
  }

  namespace internal {
    template <size_t Layer, size_t HiLayer,typename T, size_t... Ints>
    std::array<T, Layer>
    head_helper(const std::array<T, HiLayer> & arr,
                std::index_sequence<Ints...>) {
      return std::array<T, Layer> {arr[Ints]...};
    }

    template <size_t Layer, size_t HiLayer,typename T>
    std::array<T, Layer> head(const std::array<T, HiLayer> & arr) {
      return head_helper(arr, std::make_index_sequence<Layer>{});
    }

  }  // internal

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
  template<size_t Order, size_t Layer>
  class ClusterRefKey: ClusterRefBase
  {
  public:
    /**
     * Index array types need both a constant and a non-constant version. The
     * non-const version can and needs to be cast into a const version in
     * argument.
     */
    using Parent = ClusterRefBase;
    using IndexConstArray = Eigen::Map<const Eigen::Matrix<size_t, Layer+1, 1>>;
    using IndexArray = Eigen::Map<Eigen::Matrix<size_t, Layer+1, 1>>;

    //! Default constructor
    ClusterRefKey() = delete;

    /**
     * direct constructor. Initialized with an array of atoms indices,
     * and a cluster reference data
     */
    ClusterRefKey(std::array<int, Order> atom_indices,
                  IndexConstArray cluster_indices) :
      Parent{Order, Layer}, atom_indices{atom_indices},
      cluster_indices{cluster_indices.data()}
    {}

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

    const inline std::array<int, Order> & get_atom_indices() const {
      return this->atom_indices;
    }

    const int & front() const{return this->atom_indices.front();}
    const int & back() const{return this->atom_indices.back();}

    inline size_t get_cluster_index(const size_t layer) const {
      return this->cluster_indices(layer);
    }

    inline IndexConstArray get_cluster_indices() const {
      return this->cluster_indices;
    }

    constexpr static inline size_t order() {return Order;}

    constexpr static inline size_t cluster_layer() {return Layer;}

  protected:
    /**
     *  Array with unique atom indices. These can be user defined to refer to
     *  the exact same atom, e.g. in a Monte-Carlo simulation, where atoms are
     *  swapped.
     */
    std::array<int, Order> atom_indices;
    /**
     * Cluster indices by layer order, highest layer, means last
     * adaptor, and mean last entry (.back())
     */
    IndexConstArray cluster_indices;

  private:
  };

} // rascal

#endif /* CLUSTERREFBASE_H */
