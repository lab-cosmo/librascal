/**
 * file   structure_manager.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Interface for neighbourhood managers
 *
 * Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_HH_
#define SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_HH_

/**
 * Each actual implementation of a StructureManager is based on the given
 * interface
 */
#include "structure_managers/structure_manager_base.hh"
#include "structure_managers/property.hh"
#include "structure_managers/property_block_sparse.hh"
#include "structure_managers/cluster_ref_key.hh"
#include "rascal_utility.hh"
#include "json_io.hh"

//! Some data types and operations are based on the Eigen library
#include <Eigen/Dense>

//! And standard header inclusion
#include <cstddef>
#include <array>
#include <type_traits>
#include <utility>
#include <limits>
#include <tuple>
#include <sstream>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  namespace AdaptorTraits {
    //! signals if neighbours are sorted by distance
    enum class SortedByDistance : bool { yes = true, no = false };
    //! full neighbourlist or minimal neighbourlist (no permutation of clusters)
    enum class NeighbourListType { full, half };
    //! strictness of a neighbourlist with respect to a given cutoff
    enum class Strict : bool { yes = true, no = false };  //
  }  // namespace AdaptorTraits
  /* ---------------------------------------------------------------------- */

  /**
   * traits structure to avoid incomplete types in crtp Empty because it is not
   * known what it will contain
   */
  template <class ManagerImplementation>
  struct StructureManager_traits {};

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
    struct ClusterIndexPropertyComputer_Helper<
        Manager, Order, std::index_sequence<LayersHead, LayersTail...>,
        std::tuple<TupComp...>> {
      using traits = typename Manager::traits;
      constexpr static auto ActiveLayer{
          compute_cluster_layer<Order>(typename traits::LayerByOrder{})};

      using Property_t =
          Property<size_t, Order, ActiveLayer, Manager, LayersHead + 1, 1>;
      using type = typename ClusterIndexPropertyComputer_Helper<
          Manager, Order + 1, std::index_sequence<LayersTail...>,
          std::tuple<TupComp..., Property_t>>::type;
    };

    /**
     * Recursion end. Overloads helper function and is used to cycle on the
     * objects, for the last object in the list
     */
    template <typename Manager, size_t Order, size_t LayersHead,
              typename... TupComp>
    struct ClusterIndexPropertyComputer_Helper<Manager, Order,
                                               std::index_sequence<LayersHead>,
                                               std::tuple<TupComp...>> {
      using traits = typename Manager::traits;
      constexpr static auto ActiveLayer{
          compute_cluster_layer<Order>(typename traits::LayerByOrder{})};

      using Property_t =
          Property<size_t, Order, ActiveLayer, Manager, LayersHead + 1, 1>;
      using type = std::tuple<TupComp..., Property_t>;
    };

    //! Overloads the base function to call the helper function
    template <typename Manager, size_t... Layers>
    struct ClusterIndexPropertyComputer<Manager,
                                        std::index_sequence<Layers...>> {
      using type = typename ClusterIndexPropertyComputer_Helper<
          Manager, 1, std::index_sequence<Layers...>, std::tuple<>>::type;
    };

    /**
     * Empty template helper structure is created to construct cluster indices
     * tuples
     */
    template <typename Tup, typename Manager>
    struct ClusterIndexConstructor {};

    //! Overload to build the tuple
    template <typename... PropertyTypes, typename Manager>
    struct ClusterIndexConstructor<std::tuple<PropertyTypes...>, Manager> {
      static inline decltype(auto) make(Manager & manager) {
        return std::tuple<PropertyTypes...>(
            std::move(PropertyTypes(manager))...);
      }
    };
    /* Resolves the type of the cluster_indices to the corresponding Eigen::Map
     * object
     */
    // template<typename ClusterIndicesType, size_t Layer>
    // static inline IndexConstArray_t &
    // cast_cluster_indices_const(ClusterIndicesType cluster_indices) {
    //}
    // template<size_t Order, size_t ClusterLayer, size_t ParentLayer, size_t
    // NeighbourLayer>

  }  // namespace internal

  /* ---------------------------------------------------------------------- */
  /**
   * Base class interface for neighbourhood managers. The actual implementation
   * is written in the class ManagerImplementation, and the base class both
   * inherits from it and is templated by it. This allows for compile-time
   * polymorphism without runtime cost and is called a `CRTP
   * <https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_
   *
   * It inherits from StructureManagerbase because to provide a common interface
   * to the number of clusters and from `Updateable` to be able to update the
   * structure by using a vector of `Updateables`.
   *
   * @param ManagerImplementation
   * class implementation
   */
  template <class ManagerImplementation>
  class StructureManager : public StructureManagerBase {
   public:
    using StructureManager_t = StructureManager<ManagerImplementation>;
    using traits = StructureManager_traits<ManagerImplementation>;
    //! type used to represent spatial coordinates, etc
    using Vector_t = Eigen::Matrix<double, traits::Dim, 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using ClusterIndex_t = typename internal::ClusterIndexPropertyComputer<
        StructureManager, typename traits::LayerByOrder>::type;
    using ClusterConstructor_t =
        typename internal::ClusterIndexConstructor<ClusterIndex_t,
                                                   StructureManager_t>;

    //! helper type for Property creation
    template <typename T, size_t Order, Dim_t NbRow = 1, Dim_t NbCol = 1>
    using Property_t =
        Property<T, Order, get_layer(Order, typename traits::LayerByOrder{}),
                 StructureManager_t, NbRow, NbCol>;

    //! helper type for Property creation
    template <typename T, size_t Order>
    using TypedProperty_t =
        TypedProperty<T, Order,
                      get_layer(Order, typename traits::LayerByOrder{}),
                      StructureManager_t>;

    using Key_t = std::vector<int>;
    template <typename T, size_t Order>
    using BlockSparseProperty_t =
        BlockSparseProperty<T, Order,
                            get_layer(Order, typename traits::LayerByOrder{}),
                            StructureManager_t, Key_t>;

    //! type for the hyper parameter class
    using Hypers_t = json;

    //! Default constructor
    StructureManager()
        : cluster_indices_container{ClusterConstructor_t::make(*this)} {}

    //! Copy constructor
    StructureManager(const StructureManager & other) = delete;

    //! Move constructor
    StructureManager(StructureManager && other) = default;

    //! Destructor
    virtual ~StructureManager() = default;

    //! Copy assignment operator
    StructureManager & operator=(const StructureManager & other) = delete;

    //! Move assignment operator
    StructureManager & operator=(StructureManager && other) = default;

    virtual void update_self() = 0;

    // required for the construction of vectors, etc
    constexpr static int dim() { return traits::Dim; }

    /**
     * iterator over the atoms, pairs, triplets, etc in the manager. Iterators
     * like these can be used as indices for random access in atom-, pair,
     * ... -related properties.
     */
    template <size_t Order>
    class Iterator;
    using Iterator_t = Iterator<1>;
    friend Iterator_t;
    using iterator = Iterator_t;

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

    //! A proxy class which provides iteration access to atoms and ghost atoms
    class ProxyWithGhosts;

    //! A proxy class which provides iteration acces to only ghost atoms
    class ProxyOnlyGhosts;

    //! Get an iterator for a ClusterRef<1> to access pairs of an atom
    inline Iterator_t get_iterator_at(const size_t index,
                                      const size_t offset = 0) {
      return Iterator_t(*this, index, offset);
    }

    //! start of iterator
    inline Iterator_t begin() { return Iterator_t(*this, 0, 0); }
    //! end of iterator
    inline Iterator_t end() {
      return Iterator_t(*this, this->implementation().get_size(),
                        std::numeric_limits<size_t>::max());
    }

    //! Usage of iterator including ghosts; in case no ghost atoms exist, it is
    //! an iteration over all existing center atoms
    inline ProxyWithGhosts with_ghosts() {
      return ProxyWithGhosts{this->implementation()};
    }

    //! Usage of iterator for only ghosts, in case no ghosts exist, the iterator
    //! is empty
    inline ProxyOnlyGhosts only_ghosts() {
      return ProxyOnlyGhosts{this->implementation()};
    }

    //! i.e. number of atoms
    inline size_t size() const { return this->implementation().get_size(); }

    //! number of atoms including ghosts
    inline size_t size_with_ghosts() const {
      return this->implementation().get_size_with_ghosts();
    }

    template <size_t MaxOrder = traits::MaxOrder>
    inline std::enable_if_t<(MaxOrder > 1), bool>
    get_consider_ghost_neighbours() const {
      return this->implementation().get_consider_ghost_neighbours();
    }
    template <size_t MaxOrder = traits::MaxOrder>
    inline std::enable_if_t<not(MaxOrder > 1), bool>
    get_consider_ghost_neighbours() const {
      return false;
    }
    //! number of atoms, pairs, triplets in respective manager
    inline size_t nb_clusters(size_t order) const final {
      return this->implementation().get_nb_clusters(order);
    }

    //! returns position of an atom with index ``atom_tag``
    inline Vector_ref position(const int & atom_tag) {
      return this->implementation().get_position(atom_tag);
    }

    //! returns position of an atom with an AtomRef ``atom``
    inline Vector_ref position(const AtomRef & atom) {
      return this->implementation().get_position(atom);
    }

    //! returns the atom type (convention is atomic number, but nothing is
    //! imposed apart from being an integer
    inline const int & atom_type(const int & atom_tag) const {
      return this->implementation().get_atom_type(atom_tag);
    }

    //! returns the atom type (convention is atomic number, but nothing is
    //! imposed apart from being an integer
    inline int & atom_type(const int & atom_tag) {
      return this->implementation().get_atom_type(atom_tag);
    }

    /**
     * Attach a property to a StructureManager. A given calculated property is
     * only reasonable if connected with a structure. It is also connected with
     * a sanity check so that the naming of attached properties is unique. If a
     * property with the desired `name` already exists, a runtime error is
     * thrown.
     */
    void attach_property(const std::string & name,
                         std::shared_ptr<PropertyBase> property) {
      if (this->has_property(name)) {
        std::stringstream error{};
        error << "A property of name '" << name
              << "' has already been registered";
        throw std::runtime_error(error.str());
      }
      this->properties[name] = property;
    }

    /**
     * Helper function to check if a property with the specifier `name` has
     * already been attached.
     */
    inline bool has_property(const std::string & name) const {
      return not(this->properties.find(name) == this->properties.end());
    }

    template <typename Property_t>
    void create_property(const std::string & name) {
      auto property{std::make_shared<Property_t>(this->implementation())};
      this->attach_property(name, property);
    }

    template <typename T, size_t Order, Dim_t NbRow = 1, Dim_t NbCol = 1>
    void create_property(const std::string & name) {
      return create_property<Property_t<T, Order, NbRow, NbCol>>(name);
    }

    //! Accessor for an attached property with a specifier as a string
    std::shared_ptr<PropertyBase> get_property(const std::string & name) const {
      if (not this->has_property(name)) {
        std::stringstream error{};
        error << "No property of name '" << name << "' has been registered";
        throw std::runtime_error(error.str());
      }
      return this->properties.at(name);
    }

    /**
     *  Checks if the property type of user matches the actual stored
     *  property.
     */
    template <typename UserProperty_t>
    bool check_property_t(const std::string & name) const {
      auto property = get_property(name);
      try {
        UserProperty_t::check_compatibility(*property);
      } catch (const std::runtime_error & error) {
        return false;
      }
      return true;
    }

    /**
     * Throws an error if property type given from user does not match actual
     * property type.
     * TO(all) Is the try and catch need here ? it will throw in the respective
     * check_compatibility and we get the full stack with the debugger.
     */
    template <typename UserProperty_t>
    void validate_property_t(std::shared_ptr<PropertyBase> property) const {
      try {
        UserProperty_t::check_compatibility(*property);
      } catch (const std::runtime_error & error) {
        std::stringstream err_str{};
        err_str << "Incompatible UserProperty_t used : " << error.what();
        throw std::runtime_error(err_str.str());
      }
    }

    template <typename UserProperty_t>
    void validate_property_t(const std::string & name) const {
      auto property = this->get_property(name);
      this->template validate_property_t<UserProperty_t>(property);
    }

    // #BUG8486@(till) I made the function but I don't use it, because if you
    // keep a reference of an object as member variable, you have to initialize
    // it in the initialization list of the constructor. The property is not
    // created until the adaptor's update function is invoked. So I would need
    // to give the adaptor a dummy property object for initialization. I dont
    // think this is a nice solution.
    template <typename UserProperty_t>
    UserProperty_t &
    get_validated_property_ref(const std::string & name) const {
      return *this->get_validated_property<UserProperty_t>(name);
    }
    /*  Returns the typed property. Throws an error if property type given from
     *  user does not match actual property type.
     */
    template <typename UserProperty_t>
    std::shared_ptr<UserProperty_t>
    get_validated_property(const std::string & name) const {
      auto property = this->get_property(name);
      this->template validate_property_t<UserProperty_t>(property);
      return std::static_pointer_cast<UserProperty_t>(property);
    }

    void register_property(std::shared_ptr<PropertyBase> property,
                           const std::string & name) {
      this->properties[name] = property;
    }
    /**
     * Get a property of a given name. Create it if it does not exist.
     *
     * @tparam UserProperty_t full type of the property to return
     *
     * @param name name of the property to get
     *
     * @throw runtime_error if UserProperty_t is not compatible with property
     * of the given name
     */
    template <typename UserProperty_t>
    std::shared_ptr<UserProperty_t> get_property_ptr(const std::string & name) {
      if (this->has_property(name)) {
        auto property{this->get_property(name)};
        UserProperty_t::check_compatibility(*property);
        return std::static_pointer_cast<UserProperty_t>(property);
      } else {
        auto property{std::make_shared<UserProperty_t>(this->implementation())};
        this->register_property(property, name);
        return property;
      }
    }

    template <typename UserProperty_t>
    UserProperty_t & get_property_ref(const std::string & name) {
      return *this->template get_property_ptr<UserProperty_t>(name);
    }

    template <typename T, size_t Order, Dim_t NbRow = 1, Dim_t NbCol = 1>
    std::shared_ptr<Property_t<T, Order, NbRow, NbCol>>
    get_property(const std::string & name) const {
      return this
          ->template get_validated_property<Property_t<T, Order, NbRow, NbCol>>(
              name);
    }

    inline void set_updated_property_status(const bool & is_updated) {
      for (auto & element : this->properties) {
        auto & property{element.second};
        property->set_updated_status(is_updated);
      }
    }

    inline void set_updated_property_status(const std::string & name,
                                            const bool & is_updated) {
      this->properties[name]->set_updated_status(is_updated);
    }

    //! Get the full type of the structure manager
    static decltype(auto) get_name() {
      return internal::GetTypeName<ManagerImplementation>();
    }

    //! Create a new shared pointer to the object
    std::shared_ptr<ManagerImplementation> get_shared_ptr() {
      return this->implementation().shared_from_this();
    }

    //! Create a new weak pointer to the object
    std::weak_ptr<ManagerImplementation> get_weak_ptr() {
      return std::weak_ptr<ManagerImplementation>(this->get_shared_ptr());
    }

    template <size_t Order, size_t Layer>
    inline size_t
    get_cluster_size(const ClusterRefKey<Order, Layer> & cluster) const {
      return this->implementation().get_cluster_size_impl(cluster);
    }

    /**
     * Get atom_tag of index-th neighbour of this cluster, e.g. j-th
     * neighbour of atom i or k-th neighbour of pair i-j, etc.
     * Because this function is invoked with with ClusterRefKey<1, Layer> the
     * ParentLayer and NeighbourLayer have to be optional for the case Order
     * = 1.
     */
    template <size_t Order, size_t Layer>
    inline int
    get_neighbour_atom_tag(const ClusterRefKey<Order, Layer> & cluster,
                           size_t index) const {
      return this->implementation().get_neighbour_atom_tag(cluster, index);
    }

    //! get atom_tag of the index-th atom in manager
    inline int get_neighbour_atom_tag(const StructureManager & cluster,
                                      size_t & index) const {
      return this->implementation().get_neighbour_atom_tag(cluster, index);
    }

    //! Access to offsets for access of cluster-related properties
    template <size_t Order, size_t Layer>
    inline size_t
    get_offset(const ClusterRefKey<Order, Layer> & cluster) const {
      constexpr auto layer{
          StructureManager::template cluster_layer_from_order<Order>()};
      return cluster.get_cluster_index(layer);
    }

    //! Used for building cluster indices
    template <size_t Order>
    inline size_t get_offset(const std::array<size_t, Order> & counters) const {
      return this->implementation().get_offset_impl(counters);
    }

    /* Returns the neighbour's cluster_index of order 1 from an atomic index.
     */
    size_t get_atom_index(const int atom_tag) const {
      return this->implementation().get_atom_index(atom_tag);
    }

    template <size_t Order>
    constexpr static size_t cluster_layer_from_order() {
      static_assert(Order > 0, "Order is <1 this should not be");
      return get_layer(Order, typename traits::LayerByOrder{});
    }

    /**
     * When the underlying structure changes, all computations are potentially
     * invalid. This function triggers the setting of the statue variable to
     * `false` along the tree to the managers and the properties it holds.
     */
    void send_changed_structure_signal() final {
      this->set_updated_property_status(false);
      this->set_update_status(false);
      for (auto && child : this->children) {
        if (not child.expired()) {
          child.lock()->send_changed_structure_signal();
        }
      }
    }

   protected:
    /**
     * Update itself and send update signal to children nodes
     * Should only be used in the StructureManagerRoot
     */
    void update_children() final {
      if (not this->get_update_status()) {
        this->implementation().update_self();
        this->set_update_status(true);
      }
      for (auto && child : this->children) {
        if (not child.expired()) {
          child.lock()->update_children();
        }
      }
    }

    //! returns the current layer
    template <size_t Order>
    constexpr static size_t cluster_layer() {
      return compute_cluster_layer<Order>(typename traits::LayerByOrder{});
    }

    //! recursion end, not for use
    const std::array<int, 0> get_atom_tag_list() const {
      return std::array<int, 0>{};
    }

    /* Returns the cluster size in given order and layer. Because this function
     * is invoked with with ClusterRefKey<1, Layer> the ParentLayer and
     * NeighbourLayer have to be optional for the case Order = 1.
     */

    //! returns a reference to itself
    inline StructureManager & get_manager() { return *this; }

    //! necessary casting of the type
    inline ManagerImplementation & implementation() {
      return static_cast<ManagerImplementation &>(*this);
    }
    //! returns a reference for access of the implementation
    inline const ManagerImplementation & implementation() const {
      return static_cast<const ManagerImplementation &>(*this);
    }

    //! get an array with all atoms inside
    std::array<AtomRef, 0> get_atoms() const {
      return std::array<AtomRef, 0>{};
    }
    //! Starting array for builing container in iterator
    std::array<int, 0> get_atom_ids() const { return std::array<int, 0>{}; }

    //! recursion end, not for use
    inline std::array<size_t, 1> get_counters() const {
      return std::array<size_t, 1>{};
    }
    //! access to cluster_indices_container
    inline ClusterIndex_t & get_cluster_indices_container() {
      return this->cluster_indices_container;
    }

    //! access to cluster_indices_container
    inline const ClusterIndex_t & get_cluster_indices_container() const {
      return this->cluster_indices_container;
    }

    /**
     * Tuple which contains MaxOrder number of cluster_index lists for
     * reference with increasing layer depth. It is filled upon construction
     * of the neighbourhood manager via a
     * std::get<Order>(this->cluster_indices). Higher order are constructed
     * in adaptors accordingly via the lower level indices and a
     * Order-dependend index is appended to the array.
     */
    ClusterIndex_t cluster_indices_container;

    std::map<std::string, std::shared_ptr<PropertyBase>> properties{};
  };

  /* ----------------------------------------------------------------------
   */
  namespace internal {
    //! helper function that allows to append extra elements to an array It
    //! returns the given array, plus one element
    template <typename T, size_t Size, int... Indices>
    decltype(auto) append_array_helper(const std::array<T, Size> & arr, T && t,
                                       std::integer_sequence<int, Indices...>) {
      return std::array<T, Size + 1>{arr[Indices]..., std::forward<T>(t)};
    }

    //! template function allows to add an element to an array
    template <typename T, size_t Size>
    decltype(auto) append_array(const std::array<T, Size> & arr, T && t) {
      return append_array_helper(arr, std::forward<T>(t),
                                 std::make_integer_sequence<int, Size>{});
    }

    // #BUG8486@(till) changed name Counters to Container, because we handle an
    // object which can be as manager or clusterref, and this is is name
    // container in the iterator, counters is used for the the list of all
    // iterator indices
    /* ----------------------------------------------------------------------
     */
    /**
     * static branching to redirect to the correct function to get sizes,
     * offsets and neighbours. Used later by adaptors which modify or extend
     * the neighbourlist to access the correct offset.
     */
    template <bool AtMaxOrder>
    struct IncreaseHelper {
      template <class Manager_t, class Cluster_t>
      inline static size_t get_cluster_size(const Manager_t & /*manager*/,
                                            const Cluster_t & /*cluster*/) {
        throw std::runtime_error("This branch should never exist"
                                 " (cluster size).");
      }
      template <class Manager_t, class Container_t>
      inline static size_t get_offset(const Manager_t & /*manager*/,
                                      const Container_t & /*container*/) {
        throw std::runtime_error("This branch should never exist"
                                 " (offset implementation).");
      }
      template <class Manager_t, class Container_t>
      inline static int
      get_neighbour_atom_tag(const Manager_t & /*manager*/,
                             const Container_t & /*container*/,
                             size_t /*index*/) {
        throw std::runtime_error("This branch should never exist"
                                 " (cluster neigbour).");
      }
    };

    //! specialization for not at MaxOrder, these refer to the underlying
    //! manager
    template <>
    struct IncreaseHelper<false> {
      template <class Manager_t, class Cluster_t>
      inline static size_t get_cluster_size(const Manager_t & manager,
                                            const Cluster_t & cluster) {
        return manager.get_cluster_size_impl(cluster);
      }

      template <class Manager_t, class Container_t>
      inline static size_t get_offset(const Manager_t & manager,
                                      const Container_t & container) {
        return manager.get_offset_impl(container);
      }

      template <class Manager_t, class Container_t>
      inline static int get_neighbour_atom_tag(const Manager_t & manager,
                                               const Container_t & container,
                                               size_t index) {
        return manager.get_neighbour_atom_tag(container, index);
      }
    };

    template <class ParentClass, typename ClusterIndicesType, size_t Layer>
    struct ClusterIndicesConstCaster {
      using IndexConstArray_t = typename ParentClass::IndexConstArray;
      using IndexArray_t = typename ParentClass::IndexArray;

      template <typename ClusterIndicesType_ = ClusterIndicesType,
                typename std::enable_if_t<
                    std::is_same<ClusterIndicesType_, IndexConstArray_t>::value,
                    int> = 0>
      static inline IndexConstArray_t &
      cast(IndexConstArray_t & cluster_indices) {
        return cluster_indices;
      }
      template <
          typename ClusterIndicesType_ = ClusterIndicesType,
          typename std::enable_if_t<
              std::is_same<ClusterIndicesType_, IndexArray_t>::value, int> = 0>
      static inline IndexConstArray_t
      cast(const IndexArray_t & cluster_indices) {
        return IndexConstArray_t(cluster_indices.data());
      }
      template <typename ClusterIndicesType_ = ClusterIndicesType,
                typename std::enable_if_t<
                    std::is_same<ClusterIndicesType_, size_t>::value, int> = 0>
      static inline IndexConstArray_t cast(const size_t & cluster_index) {
        return IndexConstArray_t(&cluster_index);
      }
    };
  }  // namespace internal

  /* ----------------------------------------------------------------------
   */
  /**
   * Definition of the ``AtomRef`` class. It is the return type when
   * iterating over the first order of a manager.
   */
  template <class ManagerImplementation>
  class StructureManager<ManagerImplementation>::AtomRef {
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
    ~AtomRef() {}

    //! Copy assignment operator
    AtomRef & operator=(const AtomRef & other) = delete;

    //! Move assignment operator
    AtomRef & operator=(AtomRef && other) = default;

    //! return index of the atom
    inline const int & get_index() const { return this->index; }

    //! return position vector of the atom
    inline Vector_ref get_position() {
      return this->manager.position(this->index);
    }

    /**
     * return atom type (idea: corresponding atomic number, but is allowed
     * to be arbitrary as long as it is an integer)
     */
    inline const int & get_atom_type() const {
      return this->manager.atom_type(this->index);
    }
    /**
     * return atom type (idea: corresponding atomic number, but is allowed
     * to be arbitrary as long as it is an integer)
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
  };

  /* ----------------------------------------------------------------------
   */
  /**
   * Class definitionobject when iterating over the manager, then atoms,
   * then pairs, etc. in deeper Orders. This object itself is iterable again
   * up to the corresponding MaxOrder of the manager. I.e. iterating over a
   * manager provides atoms; iterating over atoms gives its pairs, etc.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  class StructureManager<ManagerImplementation>::ClusterRef
      : public ClusterRefKey<
            Order,
            ManagerImplementation::template cluster_layer_from_order<Order>()> {
   public:
    static_assert(Order > 0, "Order < 1 is not allowed.");
    using Manager_t = StructureManager<ManagerImplementation>;
    using traits = StructureManager_traits<ManagerImplementation>;
    constexpr static size_t ClusterLayer{
        ManagerImplementation::template cluster_layer_from_order<Order>()};
    using ThisParentClass = ClusterRefKey<Order, ClusterLayer>;

    using AtomRef_t = typename Manager_t::AtomRef;
    using Iterator_t = typename Manager_t::template Iterator<Order>;
    using Atoms_t = std::array<AtomRef_t, Order>;
    using iterator = typename Manager_t::template Iterator<Order + 1>;
    friend iterator;

    using IndexConstArray_t = typename ThisParentClass::IndexConstArray;
    using IndexArray_t = typename ThisParentClass::IndexArray;

    //! Default constructor
    ClusterRef() = delete;

    template <typename ClusterIndicesType>
    ClusterRef(Iterator_t & it, const std::array<int, Order> & atom_tag_list,
               ClusterIndicesType & cluster_indices)
        : ThisParentClass{atom_tag_list,
                          internal::ClusterIndicesConstCaster<
                              ThisParentClass, ClusterIndicesType,
                              ClusterLayer>::cast(cluster_indices)},
          it{it} {}

    /**
     * This is a ClusterRef of Order=1, constructed from a higher Order.
     * This function here is self referencing right now. A ClusterRefKey
     * with Order=1 is noeeded to construct it ?!
     */
    ClusterRef(ClusterRefKey<1, 0> & cluster, Manager_t & manager)
        : ClusterRefKey<1, 0>(cluster.get_atom_tag_list(),
                              cluster.get_cluster_indices()),
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
     * Returns the position of the last atom in the cluster, e.g. when
     * cluster order==1 it is the atom position, when cluster order==2 it is
     * the neighbour position, etc.
     */
    inline decltype(auto) get_position() {
      return this->get_manager().position(this->get_atom_tag());
    }

    //! returns the type of the last atom in the cluster
    inline int & get_atom_type() {
      auto && id{this->get_atom_tag()};
      return this->get_manager().atom_type(id);
    }

    //! returns the type of the last atom in the cluster
    inline const int & get_atom_type() const {
      auto && id{this->get_atom_tag()};
      return this->get_manager().atom_type(id);
    }

    /**
     * build a array of atom types from the atoms in this cluster
     */
    std::array<int, Order> get_atom_types() const;

    //! return the index of the atom/pair/etc. it is always the last one,
    //! since the other ones are accessed an Order above.
    inline int get_atom_tag() const { return this->back(); }
    //! returns a reference to the manager with the maximum layer
    inline Manager_t & get_manager() { return this->it.get_manager(); }

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
    inline size_t size() { return this->get_manager().get_cluster_size(*this); }
    //! return iterator index - this is used in cluster_indices_container as
    //! well as accessing properties
    inline size_t get_index() const { return this->it.index; }
    //! returns the clusters index (e.g. the 4-th pair of all pairs in this
    //! iteration)
    // TODO(alex) should be equal cluster_indices.back()
    // improve documentation
    inline size_t get_global_index() const {
      return this->get_manager().get_offset(*this);
    }
    //! returns the atom tags, which constitute the cluster
    const std::array<int, Order> & get_atom_tag_list() const {
      return this->atom_tag_list;
    }

    inline Iterator_t & get_iterator() { return this->it; }
    inline const Iterator_t & get_iterator() const { return this->it; }

   protected:
    //! counters for access
    inline std::array<size_t, 1> get_counters() const {
      return this->it.get_counters();
    }
    inline std::array<size_t, 1> get_offsets() const {
      return this->it.get_offsets();
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
                              std::index_sequence<Indices...> /*indices*/) {
      return std::array<int, Order>{manager.atom_type(array[Indices])...};
    }

    template <class Manager, size_t Order, size_t... Indices>
    std::array<int, Order>
    species_aggregator(const std::array<int, Order> & index_array,
                       const Manager & manager) {
      return species_aggregator_helper<Manager>(
          index_array, manager, std::make_index_sequence<Order>{});
    }
  };  // namespace internal

  template <class ManagerImplementation>
  template <size_t Order>
  std::array<int, Order>
  StructureManager<ManagerImplementation>::ClusterRef<Order>::get_atom_types()
      const {
    return internal::species_aggregator<StructureManager>(this->atom_tag_list,
                                                          this->get_manager());
  }

  /* ----------------------------------------------------------------------
   */
  /**
   * Helper functions to avoid needing dereferencing a manager in a
   * shared_ptr to loop over the centers.
   */
  template <typename T>
  auto inline begin(std::shared_ptr<T> ptr) -> typename T::iterator {
    return ptr->begin();
  }

  template <typename T>
  auto inline end(std::shared_ptr<T> ptr) -> typename T::iterator {
    return ptr->end();
  }

  /* ----------------------------------------------------------------------
   */
  /**
   * Class definition of the iterator. This is used by all clusters. It is
   * specialized for the case Order=1, when iterating over a manager and the
   * dereference are atoms, because then the container is a manager. For all
   * other cases, the container is the cluster of the Order below.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  class StructureManager<ManagerImplementation>::Iterator {
   public:
    using Manager_t = StructureManager<ManagerImplementation>;
    friend Manager_t;
    using ClusterRef_t = typename Manager_t::template ClusterRef<Order>;

    friend ClusterRef_t;
    // determine the container type
    using Container_t =
        std::conditional_t<Order == 1, Manager_t,
                           typename Manager_t::template ClusterRef<Order - 1>>;
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
    Iterator() = delete;

    //! Copy constructor
    Iterator(const Iterator & other) = default;

    //! Move constructor
    Iterator(Iterator && other) = default;

    //! Destructor
    virtual ~Iterator() = default;

    //! Copy assignment operator
    Iterator & operator=(const Iterator & other) = default;

    //! Move assignment operator
    Iterator & operator=(Iterator && other) = default;

    //! pre-increment
    inline Iterator & operator++() {
      ++this->index;
      return *this;
    }

    //! pre-decrement
    inline Iterator & operator--() {
      --this->index;
      return *this;
    }

    inline value_type operator*() {
      auto & cluster_indices_properties = std::get<Order - 1>(
          this->get_manager().get_cluster_indices_container());
      using Ref_t = typename std::remove_reference_t<decltype(
          cluster_indices_properties)>::reference;
      Ref_t cluster_indices =
          cluster_indices_properties[this->get_cluster_index()];
      return ClusterRef_t(*this, this->get_atom_tag_list(), cluster_indices);
    }

    inline const value_type operator*() const {
      const auto & cluster_indices_properties = std::get<Order - 1>(
          this->get_manager().get_cluster_indices_container());
      using Ref_t = typename std::remove_reference_t<decltype(
          cluster_indices_properties)>::const_reference;
      Ref_t cluster_indices =
          cluster_indices_properties[this->get_cluster_index()];

      const auto indices{this->get_atom_tag_list()};
      return ClusterRef_t(const_cast<Iterator &>(*this), indices,
                          cluster_indices);
    }

    //! equality
    inline bool operator==(const Iterator & other) const {
      return this->index == other.index;
    }

    //! inequality
    inline bool operator!=(const Iterator & other) const {
      return not(*this == other);
    }

    /**
     * const access to container
     */
    inline const Container_t & get_container() const { return this->container; }

   protected:
    //! constructor with container ref and starting point
    Iterator(Container_t & cont, size_t start, size_t offset)
        : container{cont}, index{start}, offset{offset} {}

    //! add atomic indices in current iteration
    std::array<int, Order> get_atom_tag_list() {
      return internal::append_array(
          container.get_atom_tag_list(),
          this->get_manager().get_neighbour_atom_tag(container, this->index));
    }

    //! add atomic indices in current iteration
    std::array<int, Order> get_atom_tag_list() const {
      return internal::append_array(
          container.get_atom_tag_list(),
          this->get_manager().get_neighbour_atom_tag(container, this->index));
    }

    //! returns the current index of the cluster in iteration
    inline size_t get_cluster_index() const {
      return this->index + this->offset;
    }

    //! returns a reference to the underlying manager at every Order
    inline Manager_t & get_manager() { return this->container.get_manager(); }
    //! returns a const reference to the underlying manager at every Order
    inline const Manager_t & get_manager() const {
      return this->container.get_manager();
    }

    //! returns the counters - which is the position in a list at each
    //! Order.
    inline std::array<size_t, Order> get_counters() {
      std::array<size_t, Order> counters;
      counters[Order - 1] = this->index;
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
    inline std::array<size_t, Order> get_offsets() {
      std::array<size_t, Order> offsets;
      offsets[Order - 1] = this->offset;
      if (Order == 1) {
        return offsets;
      } else {
        auto parental_offsets = this->container.get_offsets();
        for (size_t i{0}; i < Order - 1; i++) {
          offsets[i] = parental_offsets[i];
        }
        return offsets;
      }
    }
    //! in ascending order, this is: manager, atom, pair, triplet (i.e.
    //! cluster of Order 0, 1, 2, 3, ...
    Container_t & container;
    //! the iterators index (for moving forwards)
    size_t index;
    //! offset for access in a neighbour list during construction of the
    //! begin()
    const size_t offset;
  };

  /* ----------------------------------------------------------------------
   */
  /**
   * A class which provides the iteration range from start to the end of the
   * atoms including additional ghost atoms.
   */
  template <class ManagerImplementation>
  class StructureManager<ManagerImplementation>::ProxyWithGhosts {
   public:
    using Iterator_t =
        typename StructureManager<ManagerImplementation>::Iterator_t;

    //! Default constructor
    ProxyWithGhosts() = delete;

    //! Constructor
    ProxyWithGhosts(ManagerImplementation & manager) : manager{manager} {};

    //! Copy constructor
    ProxyWithGhosts(const ProxyWithGhosts & other) = delete;

    //! Move constructor
    ProxyWithGhosts(ProxyWithGhosts && other) = default;

    //! Destructor
    virtual ~ProxyWithGhosts() = default;

    //! Copy assignment operator
    ProxyWithGhosts & operator=(const ProxyWithGhosts & other) = delete;

    //! Move assignment operator
    ProxyWithGhosts & operator=(ProxyWithGhosts && other) = default;

    //! Start of atom list
    inline Iterator_t begin() { return this->manager.begin(); }

    //! End is all atoms including ghosts
    inline Iterator_t end() {
      return Iterator_t(this->manager, this->manager.size_with_ghosts(),
                        std::numeric_limits<size_t>::max());
    }

   protected:
    ManagerImplementation & manager;

   private:
  };

  /* ----------------------------------------------------------------------
   */
  /**
   * A class which provides the iteration range from for all ghost atoms in
   * the structure. If no ghost atoms exist, the iterator is of size zero.
   */
  template <class ManagerImplementation>
  class StructureManager<ManagerImplementation>::ProxyOnlyGhosts
      : public StructureManager<ManagerImplementation>::ProxyWithGhosts {
   public:
    using Iterator_t =
        typename StructureManager<ManagerImplementation>::Iterator_t;
    using Parent =
        typename StructureManager<ManagerImplementation>::ProxyWithGhosts;

    //! Default constructor
    ProxyOnlyGhosts() = delete;

    //!
    ProxyOnlyGhosts(ManagerImplementation & manager) : Parent{manager} {};

    //! Copy constructor
    ProxyOnlyGhosts(const ProxyOnlyGhosts & other) = delete;

    //! Move constructor
    ProxyOnlyGhosts(ProxyOnlyGhosts && other) = default;

    //! Destructor
    virtual ~ProxyOnlyGhosts() = default;

    //! Copy assignment operator
    ProxyOnlyGhosts & operator=(const ProxyOnlyGhosts & other) = delete;

    //! Move assignment operator
    ProxyOnlyGhosts & operator=(ProxyOnlyGhosts && other) = default;

    //! Start iteration at first ghost atom
    inline Iterator_t begin() { return this->manager.end(); }

   protected:
   private:
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_HH_
