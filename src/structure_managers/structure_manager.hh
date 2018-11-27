/**
 * file   structure_manager.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Interface for neighbourhood managers
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

//! header guards
#ifndef STRUCTURE_MANAGER_H
#define STRUCTURE_MANAGER_H

/**
 * Each actual implementation of a StructureManager is based on the given
 * interface
 */
#include "structure_managers/structure_manager_base.hh"
#include "structure_managers/property.hh"
#include "structure_managers/cluster_ref_key.hh"

//! Some data types and operations are based on the Eigen library
#include <Eigen/Dense>

//! And standard header inclusion
#include <cstddef>
#include <array>
#include <type_traits>
#include <utility>
#include <limits>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  namespace AdaptorTraits {
    //! signals if neighbours are sorted by distance
    enum class SortedByDistance: bool {yes = true, no = false};
    //! full neighbourlist or minimal neighbourlist (no permutation of clusters)
    enum class NeighbourListType {full, half};
    //! strictness of a neighbourlist with respect to a given cutoff
    enum class Strict:bool {yes = true, no = false}; //

    // TODO: needed?
    class Type; // type_id
  }  // AdaptorTraits
  /* ---------------------------------------------------------------------- */

  /**
   * traits structure to avoid incomplete types in crtp Empty because it is not
   * known what it will contain
   */
  template <class ManagerImplementation>
  struct StructureManager_traits
  {};

  /* ---------------------------------------------------------------------- */
  namespace internal {
    /**
     * Helper function to calculate cluster_indices_container by layer.  An
     * empty template structure is created to manage the clusters and their
     * layer
     */
    template <typename Manager, typename sequence>
    struct ClusterIndexPropertyComputer {};

    /**
     * Empty template helper structure is created to manage the clusters and
     * their layer
     */
    template <typename Manager, size_t Order, typename sequence, typename Tup>
    struct ClusterIndexPropertyComputer_Helper {};

    /**
     * Overloads helper function and is used to cycle on the objects, returning
     * the next object in line
     */
    template <typename Manager, size_t Order, size_t LayersHead,
              size_t... LayersTail, typename... TupComp>
    struct ClusterIndexPropertyComputer_Helper<Manager, Order,
                                               std::index_sequence
                                               <LayersHead, LayersTail...>,
                                               std::tuple<TupComp...>> {
      using traits = typename Manager::traits;
      constexpr static auto ActiveLayer{
        compute_cluster_layer<Order>(typename traits::LayerByOrder{})};

      using Property_t = Property<size_t, Order, ActiveLayer, LayersHead+1, 1>;
      using type = typename ClusterIndexPropertyComputer_Helper
        <Manager, Order+1, std::index_sequence<LayersTail...>,
         std::tuple<TupComp..., Property_t>>::type;
    };

    /**
     * Recursion end. Overloads helper function and is used to cycle on the
     * objects, for the last object in the list
     */
    template <typename Manager, size_t Order, size_t LayersHead,
              typename... TupComp>
    struct ClusterIndexPropertyComputer_Helper<Manager,
                                               Order,
                                               std::index_sequence<LayersHead>,
                                               std::tuple<TupComp...>> {
      using traits = typename Manager::traits;
      constexpr static auto ActiveLayer{
        compute_cluster_layer<Order>(typename traits::LayerByOrder{})};

      using Property_t = Property<size_t, Order, ActiveLayer, LayersHead+1, 1>;
      using type = std::tuple<TupComp..., Property_t>;
    };

    //! Overloads the base function to call the helper function
    template <typename Manager, size_t... Layers>
    struct ClusterIndexPropertyComputer<Manager,
                                        std::index_sequence<Layers...>> {
      using type =
        typename
        ClusterIndexPropertyComputer_Helper<Manager, 1,
                                            std::index_sequence<Layers...>,
                                            std::tuple<>>::type;
    };

    /**
     * Empty template helper structure is created to construct cluster indices
     * tuples
     */
    template <typename Tup, typename Manager>
    struct ClusterIndexConstructor {};

    //! Overload  to build the tuple
    template <typename... PropertyTypes, typename Manager>
    struct ClusterIndexConstructor<std::tuple<PropertyTypes...>, Manager> {
      static inline decltype(auto) make(Manager & manager) {
        return std::tuple<PropertyTypes...>(std::move(PropertyTypes(manager))...);
      }
    };

  }  // internal

  /* ---------------------------------------------------------------------- */
  /**
   * Base class interface for neighbourhood managers. The actual implementation
   * is written in the class ManagerImplementation, and the base class both
   * inherits from it and is templated by it. This allows for compile-time
   * polymorphism without runtime cost and is called a `CRTP
   * <https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_
   *
   * @param ManagerImplementation
   * class implementation
   */
  template <class ManagerImplementation>
  class StructureManager: public StructureManagerBase
  {
  public:
    using traits = StructureManager_traits<ManagerImplementation>;
    //! type used to represent spatial coordinates, etc
    using Vector_t = Eigen::Matrix<double, traits::Dim, 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using ClusterIndex_t = typename internal::ClusterIndexPropertyComputer
      <StructureManager, typename traits::LayerByOrder>::type;
    using ClusterConstructor_t = typename internal::ClusterIndexConstructor
      <ClusterIndex_t, StructureManager>;

    //! helper type for Property creation
    template <typename T, size_t Order, Dim_t NbRow = 1, Dim_t NbCol = 1>
    using Property_t = Property<T, Order,
                                compute_cluster_layer<Order>
                                (typename traits::LayerByOrder{}),
                                NbRow, NbCol>;

    //! helper type for Property creation
    template <typename T, size_t Order>
    using TypedProperty_t = TypedProperty<T, Order,
                                          compute_cluster_layer<Order>
                                          (typename traits::LayerByOrder{})>;

    //! Default constructor
    StructureManager() :
      cluster_indices_container{ClusterConstructor_t::make(*this)} {
    }

    //! Copy constructor
    StructureManager(const StructureManager & other) = delete;

    //! Move constructor
    StructureManager(StructureManager && other) = default;

    //! Destructor
    virtual ~StructureManager() = default;

    //! Copy assignment operator
    StructureManager
    & operator=(const StructureManager & other) = delete;

    //! Move assignment operator
    StructureManager
    & operator=(StructureManager && other)  = default;

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    /**
     * iterator over the atoms, pairs, triplets, etc in the manager. Iterators
     * like these can be used as indices for random access in atom-, pair,
     * ... -related properties.
     */
    template <size_t Order>
    class iterator;
    using Iterator_t = iterator<1>;
    friend Iterator_t;

    /**
     * return type for iterators: a light-weight atom reference, giving access
     * to an atom's position and force
     */
    class AtomRef;

    /**
     * return type for iterators: a light-weight pair, triplet, etc reference,
     * giving access to the AtomRefs of all implicated atoms
     */
    template <size_t Order>
    class ClusterRef;

    //! Get an iterator for a ClusterRef<1> to access pairs of an atom
    inline Iterator_t get_iterator_at(const size_t index,
                                      const size_t offset=0) {
      return Iterator_t(*this, index, offset);
    }

    //! start of iterator
    inline Iterator_t begin() {return Iterator_t(*this, 0, 0);}
    //! end of iterator
    inline Iterator_t end() {
      return Iterator_t(*this,
                        this->implementation().get_size(),
                        std::numeric_limits<size_t>::max());}

    //! i.e. number of atoms
    inline size_t size() const {return this->implementation().get_size();}

    //! number of atoms, pairs, triplets in respective manager
    inline size_t nb_clusters(size_t cluster_size) const override final{
      return this->implementation().get_nb_clusters(cluster_size);
    }

    //! returns position of an atom with index ``atom_index``
    inline Vector_ref position(const int & atom_index) {
      return this->implementation().get_position(atom_index);
    }

    //! returns position of an atom with an AtomRef
    inline Vector_ref position(const AtomRef & atom) {
      return this->implementation().get_position(atom);
    }

    //! returns the atom type (convention is atomic number, but nothing is
    //! imposed apart from being an integer
    inline const int& atom_type(const int & atom_index) const {
      return this->implementation().get_atom_type(atom_index);
    }

    //! returns the atom type (convention is atomic number, but nothing is
    //! imposed apart from being an integer
    inline int& atom_type(const int & atom_index) {
      return this->implementation().get_atom_type(atom_index);
    }

  protected:
    //! returns the current layer
    template <size_t Order>
    constexpr static size_t cluster_layer(){
      return compute_cluster_layer<Order>(typename traits::LayerByOrder{});
    }

    //! recursion end, not for use
    const std::array<int, 0> get_atom_indices() const {
      return std::array<int, 0>{};
    }

    //! returns the cluster size in given order and layer
    template <size_t Order, size_t Layer>
    inline size_t cluster_size(ClusterRefKey<Order, Layer> & cluster) const {
      return this->implementation().get_cluster_size(cluster);
    }

    /**
     * get atom_index of index-th neighbour of this cluster, e.g. j-th neighbour
     * of atom i or k-th neighbour of pair i-j, etc.
     */
    template <size_t Order, size_t Layer>
    inline int cluster_neighbour(ClusterRefKey<Order, Layer> & cluster,
                                 size_t index) const {
      return this->implementation().get_cluster_neighbour(cluster, index);
    }

    //! get atom_index of the index-th atom in manager
    inline int cluster_neighbour(StructureManager & cluster,
                                 size_t & index) const {
      return this->implementation().get_cluster_neighbour(cluster, index);
    }

    //! returns a reference to itself
    inline StructureManager & get_manager() {return *this;}
    //! necessary casting of the type
    inline ManagerImplementation & implementation() {
      return static_cast<ManagerImplementation&>(*this);
    }
    //! returns a reference for access of the implementation
    inline const ManagerImplementation & implementation() const {
      return static_cast<const ManagerImplementation&>(*this);
    }
    //! get an array with all atoms inside
    std::array<AtomRef, 0> get_atoms() const {
      return std::array<AtomRef, 0>{};
    }
    //! Starting array for builing container in iterator
    std::array<int, 0> get_atom_ids() const {
      return std::array<int, 0>{};
    }
    //! Access to offsets for access of cluster-related properties
    template <size_t Order, size_t CallerLayer>
    inline size_t get_offset(const ClusterRefKey<Order,
                             CallerLayer> & cluster) const {
      constexpr auto
        layer{StructureManager::template cluster_layer<Order>()};
      return cluster.get_cluster_index(layer);
    }

    //! Used for building cluster indices
    template <size_t Order>
    inline size_t get_offset(const std::array<size_t, Order> & counters) const {
      return this->implementation().get_offset_impl(counters);
    }
    //! recursion end, not for use
    inline std::array<size_t, 1> get_counters() const {
      return std::array<size_t, 1>{};
    }
    //! access to cluster_indices_container
    inline ClusterIndex_t & get_cluster_indices_container() {
      return this->cluster_indices_container;
    }

    //! access to cluster_indices_container
    inline const ClusterIndex_t & get_cluster_indices_container() const  {
      return this->cluster_indices_container;
    }

    /**
     * Tuple which contains MaxOrder number of cluster_index lists for reference
     * with increasing layer depth. It is filled upon construction of the
     * neighbourhood manager via a
     * std::get<Order>(this->cluster_indices). Higher order are constructed in
     * adaptors accordingly via the lower level indices and a Order-dependend
     * index is appended to the array.
     */
    // TODO: possible: tuple of shared pointer. adaptor_increase_maxleve makes a
    // cluster_ref_base maps onto a column of this, which referes to the
    // current cluster
    ClusterIndex_t cluster_indices_container;

  private:
  };

  /* ---------------------------------------------------------------------- */
  namespace internal {
    //! helper function that allows to append extra elements to an array It
    //! returns the given array, plus one element
    template <typename T, size_t Size, int... Indices>
    decltype(auto) append_array_helper(const std::array<T, Size> & arr, T &&  t,
                                       std::integer_sequence<int, Indices...>) {
      return std::array<T, Size+1> {arr[Indices]..., std::forward<T>(t)};
    }

    //! template function allows to add an element to an array
    template <typename T, size_t Size>
    decltype(auto) append_array (const std::array<T, Size> & arr, T &&  t) {
      return append_array_helper(arr, std::forward<T>(t),
                                 std::make_integer_sequence<int, Size>{});
    }

    /* ---------------------------------------------------------------------- */
    /**
     * static branching to redirect to the correct function to get sizes,
     * offsets and neighbours. Used later by adaptors which modify or extend the
     * neighbourlist to access the correct offset.
     */
    template<bool AtMaxOrder>
    struct IncreaseHelper {
      template<class Manager_t, class Cluster_t>
      inline static size_t get_cluster_size(const Manager_t & /*manager*/,
                                            const Cluster_t & /*cluster*/) {
        throw std::runtime_error("This branch should never exist"
                                 " (cluster size).");
      }
      template<class Manager_t, class Counters_t>
      inline static size_t get_offset_impl(const Manager_t & /*manager*/,
                                           const Counters_t & /*counters*/) {
        throw std::runtime_error("This branch should never exist"
                                 " (offset implementation).");
      }
      template<class Manager_t, class Counters_t>
      inline static
      size_t get_cluster_neighbour(const Manager_t & /*manager*/,
                                   const Counters_t & /*counters*/,
                                   size_t /*index*/) {
        throw std::runtime_error("This branch should never exist"
                                 "(cluster neigbour).");
      }
    };

    //! specialization for not at MaxOrder, these refer to the underlying
    //! manager
    template<>
    struct IncreaseHelper<false> {

      template<class Manager_t, class Cluster_t>
      inline static size_t get_cluster_size(const Manager_t & manager,
                                            const Cluster_t & cluster) {
        return manager.get_cluster_size(cluster);
      }

      template<class Manager_t, class Counters_t>
      inline static size_t get_offset_impl(const Manager_t & manager,
                                           const Counters_t & counters) {
        return manager.get_offset_impl(counters);
      }

      template<class Manager_t, class Counters_t>
      inline static size_t get_cluster_neighbour(const Manager_t & manager,
                                                 const Counters_t & counters,
                                                 size_t index) {
        return manager.get_cluster_neighbour(counters, index);
      }
    };
  }  // internal

  /* ---------------------------------------------------------------------- */
  /**
   * Definition of the ``AtomRef`` class. It is the return type when iterating
   * over the first order of a manager.
   */
  template <class ManagerImplementation>
  class StructureManager<ManagerImplementation>::AtomRef
  {
  public:
    using Vector_t = Eigen::Matrix<double, ManagerImplementation::dim(), 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using Manager_t = StructureManager<ManagerImplementation>;

    //! Default constructor
    AtomRef() = delete;

    //! constructor from iterator
    AtomRef(Manager_t & manager, const int & id)
      : manager{manager}, index{id} {}

    //! Copy constructor
    AtomRef(const AtomRef & other) = default;

    //! Move constructor
    AtomRef(AtomRef && other) = default;

    //! Destructor
    ~AtomRef(){};

    //! Copy assignment operator
    AtomRef & operator=(const AtomRef & other) = delete;

    //! Move assignment operator
    AtomRef & operator=(AtomRef && other) = default;

    //! return index of the atom
    inline const int & get_index() const {return this->index;}

    //! return position vector of the atom
    inline Vector_ref get_position() {
      return this->manager.position(this->index);
    }

    /**
     * return atom type (idea: corresponding atomic number, but is allowed to be
     * arbitrary as long as it is an integer)
     */
    inline const int & get_atom_type() const {
      return this->manager.atom_type(this->index);
    }
    /**
     * return atom type (idea: corresponding atomic number, but is allowed to be
     * arbitrary as long as it is an integer)
     */
    inline int & get_atom_type() {
      return this->manager.atom_type(this->index);
    }

  protected:
    //! reference to the underlying manager
    Manager_t & manager;
    /**
     * The meaning of `index` is manager-dependent. There are no guaranties
     * regarding contiguity. It is used internally to absolutely address
     * atom-related properties.
     */
    int index;
  private:
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Class definitionobject when iterating over the manager, then atoms, then
   * pairs, etc. in deeper Orders. This object itself is iterable again up to
   * the corresponding MaxOrder of the manager. I.e. iterating over a manager
   * provides atoms; iterating over atoms gives its pairs, etc.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  class StructureManager<ManagerImplementation>::ClusterRef :
    public ClusterRefKey<Order,
                         ManagerImplementation::template cluster_layer<Order>()>
  {
  public:
    using Manager_t = StructureManager<ManagerImplementation>;
    using traits = StructureManager_traits<ManagerImplementation>;
    constexpr static auto ClusterLayer{
      ManagerImplementation::template cluster_layer<Order>()};
    using Parent =
      ClusterRefKey<Order,
                    ClusterLayer>;
    using AtomRef_t = typename Manager_t::AtomRef;
    using Iterator_t = typename Manager_t::template iterator<Order>;
    using Atoms_t = std::array<AtomRef_t, Order>;
    using iterator = typename Manager_t::template iterator<Order+1>;
    friend iterator;

    using IndexConstArray_t = typename Parent::IndexConstArray;
    using IndexArray_t = typename Parent::IndexArray;

    //! Default constructor
    ClusterRef() = delete;

    //! ClusterRef for multiple atoms with const IndexArray
    ClusterRef(Iterator_t & it,
               const std::array<int, Order> & atom_indices,
               const IndexConstArray_t & cluster_indices) :
      Parent{atom_indices, cluster_indices}, it{it} {}

    //! ClusterRef for multiple atoms with non const IndexArray
    ClusterRef(Iterator_t & it,
               const std::array<int, Order> & atom_indices,
               IndexArray_t & cluster_indices) :
      Parent{atom_indices, IndexConstArray_t(cluster_indices.data())}, it{it} {}

    //! ClusterRef for single atom, see `cluster_index`
    ClusterRef(Iterator_t & it,
               const std::array<int, Order> & atom_indices,
               const size_t & cluster_index) :

      Parent{atom_indices, IndexConstArray_t (& cluster_index)}, it{it} {}

    /**
     * This is a ClusterRef of Order=1, constructed from a higher Order.  This
     * function here is self referencing right now. A ClusterRefKey with
     * Order=1 is noeeded to construct it ?!
     */
    template <bool FirstOrder = (Order == 1)>
    ClusterRef(std::enable_if_t<FirstOrder, ClusterRefKey<1, 0>> & cluster,
               Manager_t & manager):
      Parent(cluster.get_atom_indices(), cluster.get_cluster_indices()),
      it(manager) {}

    //! Copy constructor
    ClusterRef(const ClusterRef & other) = delete;

    //! Move constructor
    ClusterRef(ClusterRef && other) = default;

    //! Destructor
    virtual ~ClusterRef() = default;

    //! Copy assignment operator
    ClusterRef & operator=(const ClusterRef & other) = default;

    //! Move assignment operator
    ClusterRef & operator=(ClusterRef && other) = default;

    //! return atoms of the current cluster
    const std::array<AtomRef_t, Order> & get_atoms() const {
      return this->atoms;
    }

    /**
     * Returns the position of the last atom in the cluster, e.g. when cluster
     * order==1 it is the atom position, when cluster order==2 it is the
     * neighbour position, etc.
     */
    inline decltype(auto) get_position() {
      return this->get_manager().position(this->get_atom_index());
    }

    //! returns the type of the last atom in the cluster
    inline int & get_atom_type() {
      auto && id{this->get_atom_index()};
      return this->get_manager().atom_type(id);
    }

    //! returns the type of the last atom in the cluster
    inline const int & get_atom_type() const {
      auto && id{this->get_atom_index()};
      return this->get_manager().atom_type(id);
    }

    /**
     * build a array of atom types from the atoms in this cluster
     */
    std::array<int, Order> get_atom_types() const;

    //! return the index of the atom/pair/etc. it is always the last one, since
    //! the other ones are accessed an Order above.
    inline int get_atom_index() const {
      return this->back();
    }
    //! returns a reference to the manager with the maximum layer
    inline Manager_t & get_manager() {return this->it.get_manager();}

    //! return a const reference to the manager with maximum layer
    inline const Manager_t & get_manager() const {
      return this->it.get_manager();
    }
    //! start of the iteration over the cluster itself
    inline iterator begin() {
      std::array<size_t, Order> counters{this->it.get_counters()};
      auto offset = this->get_manager().get_offset(counters);
      return iterator(*this, 0, offset);
    }
    //! end of the iterations over the cluster itself
    inline iterator end() {
      return iterator(*this, this->size(), std::numeric_limits<size_t>::max());
    }
    //! returns its own size
    inline size_t size() {return this->get_manager().cluster_size(*this);}
    //! return iterator index - this is used in cluster_indices_container as
    //! well as accessing properties
    inline size_t get_index() const {
      return this->it.index;
    }
    //! returns the clusters index (e.g. the 4-th pair of all pairs in this
    //! iteration)
    inline size_t get_global_index() const {
      return this->get_manager().get_offset(*this);
    }
    //! returns the atom indices, which constitute the cluster
    const std::array<int, Order> & get_atom_indices() const {
      return this->atom_indices;
    }

    inline Iterator_t & get_iterator() {return this->it;}
    inline const Iterator_t & get_iterator() const {return this->it;}

  protected:
    //! counters for access
    inline std::array<size_t, 1> get_counters() const {
      return this->it.get_counters();
    }
    //!`atom_cluster_indices` is an initially contiguous numbering of atoms
    Iterator_t & it;
  private:
  };


  namespace internal {
    template <class Manager, size_t Order, size_t... Indices>
    std::array<int, Order>
    species_aggregator_helper(const std::array<int, Order> & array,
                              const Manager & manager,
                              std::index_sequence<Indices...> /*indices*/){
      return std::array<int, Order>{
        manager.atom_type(array[Indices])...};
    }

    template <class Manager, size_t Order, size_t... Indices>
    std::array<int, Order>
    species_aggregator(const std::array<int, Order> &index_array,
                       const Manager & manager) {
      return species_aggregator_helper<Manager>
        (index_array,
         manager,
         std::make_index_sequence<Order>{});
    }
  }; // internal

  template <class ManagerImplementation>
  template <size_t Order>
  std::array<int, Order> StructureManager<ManagerImplementation>::
  ClusterRef<Order>::get_atom_types() const
   {
     return internal::species_aggregator<StructureManager>(this->atom_indices,
                                                           this->get_manager());
   }

  /* ---------------------------------------------------------------------- */
  /**
   * Class definition of the iterator. This is used by all clusters. It is
   * specialized for the case Order=1, when iterating over a manager and the
   * dereference are atoms, because then the container is a manager. For all
   * other cases, the container is the cluster of the Order below.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  class StructureManager<ManagerImplementation>::iterator
  {
  public:
    using Manager_t = StructureManager<ManagerImplementation>;
    friend Manager_t;
    using ClusterRef_t = typename Manager_t::template ClusterRef<Order>;

    friend ClusterRef_t;
    // determine the container type
    using Container_t =
      std::conditional_t
      <Order == 1,
       Manager_t,
       typename Manager_t::template
       ClusterRef<Order-1>>;
    static_assert(Order > 0, "Order has to be positive");

    static_assert(Order <= traits::MaxOrder,
                  "Order > MaxOrder, impossible iterator");

    using AtomRef_t = typename Manager_t::AtomRef;

    using value_type = ClusterRef_t;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    // using iterator_category = std::random_access_iterator_tag;

    using reference = value_type;

    //! Default constructor
    iterator() = delete;

    //! Copy constructor
    iterator(const iterator & other) = default;

    //! Move constructor
    iterator(iterator && other) = default;

    //! Destructor
    virtual ~iterator() = default;

    //! Copy assignment operator
    iterator & operator=(const iterator & other) = default;

    //! Move assignment operator
    iterator & operator=(iterator && other) = default;

    //! pre-increment
    inline iterator & operator ++ () {
      ++this->index;
      return *this;
    }

    //! pre-decrement
    inline iterator & operator -- () {
      --this->index;
      return *this;
    }

    //! dereference: calculate cluster indices
    inline value_type operator * () {
      auto & cluster_indices_properties =
        std::get<Order-1>(this->get_manager().get_cluster_indices_container());
      using Ref_t = typename
        std::remove_reference_t<decltype(cluster_indices_properties)>::reference;
      Ref_t cluster_indices =
        cluster_indices_properties[this->get_cluster_index()];
      return ClusterRef_t(*this, this->get_atom_indices(), cluster_indices);
    }

    //! dereference: calculate cluster indices
    inline const value_type operator * () const {
      const auto & cluster_indices_properties =
        std::get<Order-1>(this->get_manager().get_cluster_indices_container());
      using Ref_t = typename
        std::remove_reference_t<decltype(cluster_indices_properties)>::
        const_reference;
      Ref_t cluster_indices =
        cluster_indices_properties[this->get_cluster_index()];
      const auto indices{this->get_atom_indices()};
      return ClusterRef_t(const_cast<iterator&>(*this), indices, cluster_indices);
    }

    //! equality
    inline bool operator == (const iterator & other) const {
      return this->index == other.index;
    }

    //! inequality
    inline bool operator != (const iterator & other) const {
      return not (*this == other);
    }

    /**
     * const access to container
     */
    inline const Container_t & get_container() const {return this->container;}

  protected:
    //! constructor with container ref and starting point
    iterator(Container_t & cont, size_t start, size_t offset)
      :container{cont}, index{start}, offset{offset} {}

    //! add atomic indices in current iteration
    std::array<int, Order> get_atom_indices() {
      return internal::append_array
        (container.get_atom_indices(),
         this->get_manager().cluster_neighbour(container, this->index));
    }

    //! add atomic indices in current iteration
    std::array<int, Order> get_atom_indices() const {
      return internal::append_array
        (container.get_atom_indices(),
         this->get_manager().cluster_neighbour(container, this->index));
    }

    //! returns the current index of the cluster in iteration
    inline size_t get_cluster_index() const {
      return this->index + this->offset;
    }
    //! returns a reference to the underlying manager at every Order
    inline Manager_t & get_manager() {return this->container.get_manager();}
    //! returns a const reference to the underlying manager at every Order
    inline const Manager_t & get_manager() const {
      return this->container.get_manager();
    }

    //! returns the counters - which is the position in a list at each Order.
    inline std::array<size_t, Order> get_counters() {
      std::array<size_t, Order> counters;
      counters[Order-1] = this->index;
      if (Order == 1) {
        return counters;
      } else {
        auto parental_counters = this->container.get_counters();
        for (size_t i{0}; i < Order - 1; i++) {
          counters[i] = parental_counters[i];
        }
        return counters;
      }
    }
    //! in ascending order, this is: manager, atom, pair, triplet (i.e. cluster
    //! of Order 0, 1, 2, 3, ...
    Container_t & container;
    //! the iterators index (for moving forwards)
    size_t index;
    //! offset for access in a neighbour list during construction of the begin()
    const size_t offset;
  private:
  };
}  // rascal

#endif /* STRUCTURE_MANAGER_H */
