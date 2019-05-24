/**
 * file   cluster_ref_key.hh
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

#include "structure_managers/cluster_ref_cluster_key.hh"

#include <Eigen/Dense>

#include <tuple>
#include <array>
#include <iostream>

namespace rascal {
  /* ---------------------------------------------------------------------- */
  /**
   * Layer calculations and manipulations
   */
  //! Computes layer by cluster dimension for new adaptor layer, depending on
  //! existing layer by order.
  template <size_t MaxOrder, class T>
  struct LayerIncreaser {};

  template <size_t MaxOrder, size_t... Ints>
  struct LayerIncreaser<MaxOrder, std::index_sequence<Ints...>> {
    using type = std::index_sequence<(Ints + 1)...>;
  };

  template <size_t MaxOrder, size_t... Ints>
  using LayerIncreaser_t =
      typename LayerIncreaser<MaxOrder, std::index_sequence<Ints...>>::type;

  /* ---------------------------------------------------------------------- */
  /**
   * Dynamic access to depth by cluster order
   */

  template <size_t MaxLevel, size_t... Ints>
  constexpr size_t get_depth(size_t index, std::index_sequence<Ints...>) {
    constexpr size_t arr[]{Ints...};
    return arr[index];
  }

  /* ---------------------------------------------------------------------- */
  //! Extends layer by cluster for an additional cluster dimension
  template <size_t MaxOrder, class T>
  struct LayerExtender {};

  template <size_t MaxOrder, size_t... Ints>
  struct LayerExtender<MaxOrder, std::index_sequence<Ints...>> {
    using type = std::index_sequence<Ints..., 0>;
  };

  template <size_t MaxOrder, size_t... Ints>
  using LayerExtender_t =
      typename LayerExtender<MaxOrder, std::index_sequence<Ints...>>::type;

  /* ---------------------------------------------------------------------- */
  //! Dynamic access to all layers by cluster dimension (probably not necessary)
  template <size_t MaxOrder, size_t... Ints>
  constexpr std::array<size_t, MaxOrder>
  get_layers(std::index_sequence<Ints...>) {
    return std::array<size_t, MaxOrder>{Ints...};
  }

  /* ---------------------------------------------------------------------- */
  //! extractors helpers for cluster layers
  namespace internal {
    constexpr static int InvalidLayer = -1;
    template <size_t head, size_t... tail>
    struct Min {
      constexpr static size_t value{
          head < Min<tail...>::value ? head : Min<tail...>::value};
    };

    template <size_t head>
    struct Min<head> {
      constexpr static size_t value{head};
    };

    template <class Sequence>
    struct MinExtractor {};

    template <size_t... Ints>
    struct MinExtractor<std::index_sequence<Ints...>> {
      constexpr static size_t value{Min<Ints...>::value};
    };

    template <size_t Order, class Sequence, size_t... Ints>
    struct HeadExtractor {};

    template <size_t... seq>
    struct HeadExtractorTail {
      using type = std::index_sequence<seq...>;
    };

    template <size_t Order, size_t head, size_t... tail, size_t... seq>
    struct HeadExtractor<Order, std::index_sequence<seq...>, head, tail...> {
      using Extractor_t = std::conditional_t<
          (Order > 1),
          HeadExtractor<Order - 1, std::index_sequence<seq..., head>, tail...>,
          HeadExtractorTail<seq..., head>>;
      using type = typename Extractor_t::type;
    };
  }  // namespace internal

  /* ---------------------------------------------------------------------- */
  //! returns the cluster layer for accessing properties at a specific layer in
  //! a stack
  template <size_t Order, size_t... Ints>
  constexpr size_t compute_cluster_layer(const std::index_sequence<Ints...> &) {
    using ActiveDimensions =
        typename internal::HeadExtractor<Order, std::index_sequence<>,
                                         Ints...>::type;
    return internal::MinExtractor<ActiveDimensions>::value;
  }

  // TODO(alex) check if I changed this correctly
  template <size_t... Ints>
  constexpr size_t get_layer(const size_t order,
      const std::index_sequence<Ints...>) {
    constexpr size_t arr[]{Ints...};
    return arr[order-1];
  }
  //! Dynamic access to layer by cluster dimension (possibly not necessary)
  //template <size_t MaxOrder, size_t... Ints>
  //constexpr size_t get_layer(size_t index, std::index_sequence<Ints...>) {
  //  constexpr size_t arr[]{Ints...};
  //  return arr[index];
  //}

  /**
   * Static access to layer by cluster dimension (e.g., for defining template
   * parameter `NbRow` for a property
   */
  template <size_t index, size_t... Ints>
  constexpr size_t get(std::index_sequence<Ints...>) {
    return get<index>(std::make_tuple(Ints...));
  }

  /* ---------------------------------------------------------------------- */
  namespace internal {
    //! extracts the head of the layer by order
    template <size_t Layer, size_t HiLayer, typename T, size_t... Ints>
    std::array<T, Layer> head_helper(const std::array<T, HiLayer> & arr,
                                     std::index_sequence<Ints...>) {
      return std::array<T, Layer>{arr[Ints]...};
    }
    //! specialization of the head extractor
    template <size_t Layer, size_t HiLayer, typename T>
    std::array<T, Layer> head(const std::array<T, HiLayer> & arr) {
      return head_helper(arr, std::make_index_sequence<Layer>{});
    }
  }  // namespace internal

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

  
  template <size_t Order_, size_t Layer_, size_t ParentLayer_, size_t NeighbourLayer_>
  struct ClusterRefKeyInfo {
    static_assert(not(Order_==1) || Layer_==ParentLayer_ , "Layer does not argee with ParentLayer for Order 1 ClusterRefKeyInfo");
    static_assert(not(Order_==1) || Layer_==NeighbourLayer_, "Layer does not argee with NeighbourLayer for Order 1 ClusterRefKeyInfo");
    constexpr static size_t Order = Order_; 
    constexpr static size_t Layer = Layer_; 
    constexpr static size_t ParentLayer = ParentLayer_; 
    constexpr static size_t NeighbourLayer = NeighbourLayer_; 
  };

  // ValidateDefaultTemplateParamater
  template <size_t Order, size_t Layer>
  struct ClusterRefKeyDefaultTemplateParamater {
    static_assert(Order==1, "ClusterRefKey of Order > 1 requires ParentLayer" 
        "and NeigbhourLayer. The usage ClusterRefKey<Order, Layer> is only valid"
        "for Order = 1"); 
    constexpr static size_t ParentLayer = Layer;
    constexpr static size_t NeighbourLayer = Layer;
  };

  // TODO(till) 
  /* Usage of ClusterRefKey<Order,Layer> with Order>1 is invalid usage. ClusterRefKey<1,Layer> is valid usage. ClusterRefKey<1,Layer,y,z> is valid usage, but y and z are ignored.
   */
  template <size_t Order, size_t Layer,
      size_t ParentLayer=
          ClusterRefKeyDefaultTemplateParamater<Order, Layer>::ParentLayer,
      size_t NeighbourLayer=
          ClusterRefKeyDefaultTemplateParamater<Order, Layer>::NeighbourLayer>
  class ClusterRefKey : public ClusterRefClusterKey<Order, Layer> {
   public:
    /**
     * Index array types need both a constant and a non-constant version. The
     * non-const version can and needs to be cast into a const version in
     * argument.
     */

    static_assert(not(Order==1) || Layer==ParentLayer , "Layer does not argee with ParentLayer for Order 1 ClusterRefKeyInfo");
    static_assert(not(Order==1) || Layer==NeighbourLayer, "Layer does not argee with NeighbourLayer for Order 1 ClusterRefKeyInfo");

    using ParentClass = ClusterRefClusterKey<Order, Layer>;
    using IndexConstArray = typename ParentClass::IndexConstArray;
    using IndexArray = typename ParentClass::IndexArray;
    using AtomIndex_t = typename ParentClass::AtomIndex_t;
    using NeighbourIndexConstArray =
        Eigen::Map<const Eigen::Matrix<size_t, NeighbourLayer + 1, 1>>;


    /* For Order > 1 this object contains the neighbour object, because the neighbour object is nowhere else stored, but for Order = 1 it is a reference to itself, because an object cannot store itself.
     */
    typedef ParentClass & Neighbour_t;
    typedef ParentClass Parent_t;
    /* #ParentFeature #NeighbourFeature
     *using NeighbourClusterRefKey_t = ClusterRefClusterKey<1, NeighbourLayer>;
     *typedef typename std::conditional<Order==1,
     *         ParentClass &, NeighbourClusterRefKey_t
     *             >::type Neighbour_t;
     *using ParentClusterRefKey_t =
     *    ClusterRefClusterKey<(Order==1) ? 1 : Order-1, ParentLayer>;
     *typedef typename std::conditional<Order==1,
     *         ParentClass, ParentClusterRefKey_t
     *             >::type Parent_t;
     */
    //! Default constructor
    ClusterRefKey() = delete;

    /* #ParentFeature #NeighbourFeature
     * Definition for ClusterRefKey using neighbours and parents for order > 1,
     * if we will extend this feature the constructor might be required.
     *template<size_t Order_=Order, typename std::enable_if_t<not(Order_==1),int> = 0>
     *ClusterRefKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices,
     *    Parent_t & cluster_parent,
     *    NeighbourIndexConstArray cluster_neighbour_cluster_indices) 
     *    : ParentClass{atom_indices, cluster_indices},
     *      cluster_parent{cluster_parent},
     *      cluster_neighbour{
     *          Neighbour_t(
     *              std::array<int, 1>{atom_indices.back()},
     *              // TODO(alex) in ClusterRef we do different casting for the indices
     *              // here it seems to work without any difference between
     *              // cases, this seems fishy
     *              NeighbourIndexConstArray(cluster_neighbour_cluster_indices)) 
     *      }
     *      {}
     */

    // If we extend parent and neighbour feature, this constructor is used for Order 1.
    // Currently it is used for every case.
    //template<size_t Order_=Order, typename std::enable_if_t<(Order_==1),int> = 0>
    ClusterRefKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices)
        : ParentClass{atom_indices, cluster_indices},
          cluster_parent{*this},
            cluster_neighbour{*this} {}

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

    //! returns the cluster's index, given a specific layer
    size_t get_internal_neighbour_cluster_index(const size_t layer) {
      return this->cluster_neighbour.get_cluster_index(layer);
    }
    NeighbourIndexConstArray get_internal_neighbour_cluster_indices() const {
      return this->cluster_neighbour.get_cluster_indices();
    }

    //TODO(alex) rename by using the above
    template<size_t Order_=Order>
    std::enable_if_t<not(Order_==1), size_t> get_neighbour_cluster_index(const size_t layer) const {
      return this->cluster_neighbour.get_cluster_index(layer);
    }

    template<size_t Order_=Order>
    std::enable_if_t<(Order_==1), size_t> get_neighbour_cluster_index(const size_t layer) const {
      return this->get_cluster_index(layer);
      // TODO(alex) neighbour does not seem do be pointer type which it should be
      //return this->cluster_neighbour->get_cluster_index(layer);
    }
    template<size_t Order_=Order>
    std::enable_if_t<not(Order_==1), NeighbourIndexConstArray> get_neighbour_cluster_indices() const {
      return this->cluster_neighbour.get_cluster_indices();
    }

    template<size_t Order_=Order>
    std::enable_if_t<Order_==1, IndexConstArray> get_neighbour_cluster_indices() const {
      return this->cluster_indices;
    }

   protected:
    /** If the cluster represents a tuple of form (i,j,...,k,l) the corresponding parent cluster represents the tuple (i,j,...,k) e.g for a triplet (i,j,k) we save the reference to the pair (i,j). To access information of cluster<Order-1>, because the cluster indices of a cluster<Order> are independent from that of the corresponding cluster<Order-1>.
     */
    Parent_t & cluster_parent;
    /* The internal neighbour of the cluster. The atom with the last atom index in atom_indices. For Order = 1 this will be a reference to itself, and for Order > 1 it stores the neighbour object.
     */
    Neighbour_t cluster_neighbour;
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_CLUSTER_REF_KEY_HH_
