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

#include "structure_managers/cluster_ref_base.hh"

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
    constexpr static size_t Order = Order_; 
    constexpr static size_t Layer = Layer_; 
    constexpr static size_t ParentLayer = ParentLayer_; 
    constexpr static size_t NeighbourLayer = NeighbourLayer_; 
  };

  // ValidateDefaultTemplateParamater
  template <size_t Order, size_t Layer>
  struct ClusterRefKeyDefaultTemplateParamater {
    static_assert(Order==1, "ClusterRefKey of Order > 1 requires ParentLayer and NeigbhourLayer. The usage ClusterRefKey<Order, Layer> is only valid for Order = 1"); 
    using ParentInfo_t = ClusterRefKeyInfo<1, Layer, Layer, Layer>;
    // TODO(alex) diff in usage?
    //typedef ClusterRefKeyInfo<Order, Layer, Layer, Layer> ParentInfo_t;
    constexpr static size_t NeighbourLayer = Layer;
  };

  // TODO(alex) change validate_crk_template_parameters in SM file
  // TODO(till) 
  /* Usage of ClusterRefKey<Order,Layer> with Order>1 is invalid usage. ClusterRefKey<1,Layer> is valid usage. ClusterRefKey<1,Layer,y,z> is valid usage, but y and z are ignored.
   */
  template <size_t Order, size_t Layer,
           typename ParentInfo_t = typename
               ClusterRefKeyDefaultTemplateParamater<Order, Layer>::ParentInfo_t,
           size_t NeighbourLayer =
               ClusterRefKeyDefaultTemplateParamater<Order, Layer>::NeighbourLayer>
  class ClusterRefKey : public ClusterRefBase {
   public:
    /**
     * Index array types need both a constant and a non-constant version. The
     * non-const version can and needs to be cast into a const version in
     * argument.
     */
    constexpr static size_t ParentLayer = ParentInfo_t::Layer;
    constexpr static size_t ParentNeighbourLayer = ParentInfo_t::NeighbourLayer;
    constexpr static size_t ParentParentLayer = ParentInfo_t::ParentLayer;

    using Parent = ClusterRefBase;
    using IndexConstArray =
        Eigen::Map<const Eigen::Matrix<size_t, Layer + 1, 1>>;
    using IndexArray = Eigen::Map<Eigen::Matrix<size_t, Layer + 1, 1>>;
    // TODO(alex check if the two neighbour are still needed
    using NeighbourIndexConstArray =
        Eigen::Map<const Eigen::Matrix<size_t, NeighbourLayer + 1, 1>>;
    using NeighbourIndexArray = Eigen::Map<Eigen::Matrix<
        size_t, NeighbourLayer + 1, 1>>;
    using AtomIndex_t = std::array<int, Order>;
    /* For Order > 1 this object contains the neighbour object, because the neighbour object is nowhere else stored, but for Order = 1 it is a reference to itself, because an object cannot store itself.
     */
    using ThisClusterRefKey_t = ClusterRefKey<Order, Layer, ParentInfo_t, NeighbourLayer>;
    using NeighbourClusterRefKey_t = ClusterRefKey<1, NeighbourLayer>;
    typedef typename std::conditional<Order==1,
             ThisClusterRefKey_t &,
             NeighbourClusterRefKey_t
                 >::type Neighbour_t;
    //typedef ClusterRefKeyInfo<
    //    (Order==1) ? 1 : Order-1, ParentLayer,
    //    ParentParentLayer, ParentNeighbourLayer> ParentParentInfo_t;
    using ParentParentInfo_t = ClusterRefKeyInfo<
        (Order<3) ? 1 : Order-2, ParentLayer,
        ParentParentLayer, ParentNeighbourLayer>;
    using ParentClusterRefKey_t = ClusterRefKey<(Order==1) ? 1 : Order-1,
          ParentLayer, ParentParentInfo_t, ParentNeighbourLayer>;
    typedef typename std::conditional<Order==1,
             ThisClusterRefKey_t,
             ParentClusterRefKey_t
                 >::type Parent_t;
    // Order-1 must be resolved and returns error
    //typedef typename std::conditional<Order==1, ClusterRefKey<1, Layer>, ClusterRefKey<Order-1, Layer, ParentLayer, NeighbourLayer>>::type Parent_t;
    //! Default constructor
    ClusterRefKey() = delete;

    /* Constructor for Order > 1. To use enable_if an template parameter has to be given, therefore the Order is copied. The template arguments of the construction are not in accesible by the user, therefore it cannot be invalidly used.
     */
    // seems ill defined because object is deleted on creation
    // template<size_t Order_=Order, typename std::enable_if_t<not(Order_==1),int> = 0>
    // ClusterRefKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices,
    //     Parent_t & cluster_parent,
    //     Neighbour_t cluster_neighbour) 
    //     : Parent{Order, Layer}, atom_indices{atom_indices}, 
    //       cluster_indices{cluster_indices.data()},
    //       cluster_parent{cluster_parent}, cluster_neighbour{cluster_neighbour} {}

    template<size_t Order_=Order, typename std::enable_if_t<not(Order_==1),int> = 0>
    ClusterRefKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices,
        Parent_t & cluster_parent,
        NeighbourIndexConstArray cluster_neighbour_cluster_indices) 
        : Parent{Order, Layer}, atom_indices{atom_indices}, 
          cluster_indices{cluster_indices.data()},
          cluster_parent{cluster_parent},
          cluster_neighbour{ 
              ClusterRefKey<1, NeighbourLayer>(
                  std::array<int, 1>{atom_indices.back()},
                  cluster_neighbour_cluster_indices)
          } 
          {}

    /* Constructor for Order 1 
     */
    template<size_t Order_=Order, typename std::enable_if_t<(Order_==1),int> = 0>
    ClusterRefKey(AtomIndex_t atom_indices, IndexConstArray cluster_indices)
        : Parent{Order, Layer}, atom_indices{atom_indices},           
          cluster_indices{cluster_indices.data()},
          cluster_parent{*this}, cluster_neighbour{*this} {}

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

    //! returns the atom indices of the current cluster
    const inline AtomIndex_t & get_atom_indices() const {
      return this->atom_indices;
    }

    //! returns the first atom index in this cluster
    const int & front() const { return this->atom_indices.front(); }
    //! returns the last atom index in this cluster
    const int & back() const { return this->atom_indices.back(); }

    //! returns the cluster's index, given a specific layer
    inline size_t get_cluster_index(const size_t layer) const {
      return this->cluster_indices(layer);
    }

    inline size_t get_neighbour_cluster_index(const size_t layer) const {
      return this->cluster_neighbour.get_cluster_index(layer);
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

    /** If the cluster represents a tuple of form (i,j,...,k,l) the corresponding parent cluster represents the tuple (i,j,...,k) e.g for a triplet (i,j,k) we save the reference to the pair (i,j). To access information of cluster<Order-1>, because the cluster indices of a cluster<Order> are independent from that of the corresponding cluster<Order-1>.
     */
    Parent_t & cluster_parent;
    /* The internal neighbour of the cluster. The atom with the last atom index in atom_indices. For Order = 1 this will be a reference to itself, and for Order > 1 it stores the neighbour object.
     */
    Neighbour_t cluster_neighbour;
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_CLUSTER_REF_KEY_HH_
