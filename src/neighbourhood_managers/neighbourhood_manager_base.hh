/**
 * file   neighbourhood_manager_base.hh
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

#ifndef NEIGHBOURHOOD_MANAGER_BASE_H
#define NEIGHBOURHOOD_MANAGER_BASE_H

#include "neighbourhood_managers/cluster_ref_base.hh"
#include "neighbourhood_managers/property.hh"

#include <Eigen/Dense>

#include <cstddef>
#include <array>
#include <type_traits>
#include <utility>
#include <limits>

namespace rascal {



  namespace AdaptorTraits {

    enum class SortedByDistance: bool {yes = true, no = false};
    enum class MinImageConvention: bool {yes = true, no = false};
    enum class NeighbourListType {full, half};
    //------------------------------------------------------------------------//
    enum class Strict:bool {yes = true, no = false}; // r_cut

    class Type; // type_id
  }  // AdaptorTraits

  //! traits structure to avoid incomplete types in crtp
  template <class Manager>
  struct NeighbourhoodManager_traits
  {};

  namespace internal {
    /**
     * Helper function to calculate cluster_indices by depth.
     */
    template <typename Manager, typename sequence>
    struct ClusterIndexPropertyComputer {};

    template <typename Manager, size_t Level, typename sequence, typename Tup>
    struct ClusterIndexPropertyComputer_Helper {};

    template <typename Manager, size_t Level, size_t DepthsHead,
             size_t... DepthsTail, typename... TupComp>
    struct ClusterIndexPropertyComputer_Helper<Manager,
                                               Level,
                                               std::index_sequence
                                               <DepthsHead, DepthsTail...>,
                                               std::tuple<TupComp...>> {
      using Property_t = Property<Manager, size_t, Level, DepthsHead+1, 1>;
      using type = typename ClusterIndexPropertyComputer_Helper
        <Manager, Level+1, std::index_sequence<DepthsTail...>,
         std::tuple<TupComp..., Property_t>>::type;
    };

    // Recursion end
    template <typename Manager, size_t Level, size_t DepthsHead,
             typename... TupComp>
    struct ClusterIndexPropertyComputer_Helper<Manager,
                                               Level,
                                               std::index_sequence<DepthsHead>,
                                               std::tuple<TupComp...>> {
      using Property_t = Property<Manager, size_t, Level, DepthsHead+1, 1>;
      using type = std::tuple<TupComp..., Property_t>;
    };

    template <typename Manager, size_t... Depths>
    struct ClusterIndexPropertyComputer<Manager,
                                        std::index_sequence<Depths...>> {
      using type =
        typename
        ClusterIndexPropertyComputer_Helper<Manager, 1,
                                            std::index_sequence<Depths...>,
                                            std::tuple<>>::type;
    };

    template <typename Tup, typename Manager>
    struct ClusterIndexConstructor {};

    template <typename... PropertyTypes, typename Manager>
    struct ClusterIndexConstructor<std::tuple<PropertyTypes...>, Manager> {
      static inline decltype(auto) make(Manager & manager) {
        return std::tuple<PropertyTypes...>(std::move(PropertyTypes(manager))...);
      }
    };

  }  // internal

  /**
   * Base class interface for neighbourhood managers. The actual
   * implementation is written in the class ManagerImplementation, and
   * the base class both inherits from it and is templated by it. This
   * allows for compile-time polymorphism without runtime cost and is
   * called a `CRTP
   * <https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_
   *
   * @param ManagerImplementation
   * class implementation
   */
  template <class ManagerImplementation>
  class NeighbourhoodManagerBase
  {
  public:
    using traits = NeighbourhoodManager_traits<ManagerImplementation>;
    using Vector_t = Eigen::Matrix<double, traits::Dim, 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using ClusterIndex_t = typename internal::ClusterIndexPropertyComputer
      <NeighbourhoodManagerBase, typename traits::DepthByDimension>::type;
    using ClusterConstructor_t = typename internal::ClusterIndexConstructor
      <ClusterIndex_t, NeighbourhoodManagerBase>;

    //! Default constructor
    NeighbourhoodManagerBase() :
      cluster_indices{ClusterConstructor_t::make(*this)} {

    }

    //! Copy constructor
    NeighbourhoodManagerBase(const NeighbourhoodManagerBase & other) = delete;

    //! Move constructor
    NeighbourhoodManagerBase(NeighbourhoodManagerBase && other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerBase() = default;

    //! Copy assignment operator
    NeighbourhoodManagerBase
    & operator=(const NeighbourhoodManagerBase & other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerBase
    & operator=(NeighbourhoodManagerBase && other)  = default;

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}

    /**
     * iterator over the atoms, pairs, triplets, etc in the
     * manager. Iterators like these can be used as indices for random
     * access in atom-, pair, ... -related properties.
     */
    template <size_t Level>
    class iterator;
    using Iterator_t = iterator<1>;
    friend Iterator_t;

    /**
     * return type for iterators: a light-weight atom reference,
     * giving access to an atom's position and force
     */
    class AtomRef;

    /**
     * return type for iterators: a light-weight pair, triplet, etc reference,
     * giving access to the AtomRefs of all implicated atoms
     */
    template <size_t Level>
    class ClusterRef;

    inline Iterator_t get_iterator_at(const size_t index, const size_t offset) {
      return Iterator_t(*this, index, offset);
    }
    
    inline Iterator_t begin() {return Iterator_t(*this, 0, 0);}
    inline Iterator_t end() {
      return Iterator_t(*this,
                        this->implementation().size(),
                        std::numeric_limits<size_t>::max());}
    inline size_t size() const {return this->implementation().get_size();}

    inline size_t nb_clusters(size_t cluster_size) const {
      return this->implementation().get_nb_clusters(cluster_size);
    }

    inline Vector_ref position(const int & atom_index) {
      return this->implementation().get_position(atom_index);
    }

    inline Vector_ref position(const AtomRef & atom) {
      return this->implementation().get_position(atom);
    }

    template <size_t Level, size_t Depth>
    inline Vector_ref
    neighbour_position(ClusterRefBase<Level, Depth> & cluster) {
      return this->implementation().get_neighbour_position(cluster);
    }

    inline int atom_type(const int & atom_index) {
      return this->implementation().get_atom_type(atom_index);
    }

  protected:
    template <size_t Level>
    constexpr static size_t cluster_depth(){
      return compute_cluster_depth<Level>(typename traits::DepthByDimension{});
    }

    //! recursion end, not for use
    const std::array<int, 0> get_atom_indices() const {
      return std::array<int, 0>{};
    }

    template <size_t Level, size_t Depth>
    inline size_t cluster_size(ClusterRefBase<Level, Depth> & cluster) const {
      return this->implementation().get_cluster_size(cluster);
    }

    inline size_t cluster_size(const int & atom_index) const {
      return this->implementation().get_cluster_size(atom_index);
    }

    //! get atom_index of index-th neighbour of this cluster, e.g. j-th
    // neighbour of atom i or k-th neighbour of pair i-j, etc.
    template <size_t Level, size_t Depth>
    inline int cluster_neighbour(ClusterRefBase<Level, Depth> & cluster,
                                 size_t index) const {
      return this->implementation().get_cluster_neighbour(cluster, index);
    }

    //! get atom_index of the index-th atom in manager
    inline int cluster_neighbour(NeighbourhoodManagerBase & cluster,
                                 size_t & index) const {
      return this->implementation().get_cluster_neighbour(cluster, index);
    }

    inline NeighbourhoodManagerBase & get_manager() {return *this;}

    inline ManagerImplementation & implementation() {
      return static_cast<ManagerImplementation&>(*this);
    }
    inline const ManagerImplementation & implementation() const {
      return static_cast<const ManagerImplementation&>(*this);
    }

    std::array<AtomRef, 0> get_atoms() const {
      return std::array<AtomRef, 0>{};
    }

    //! Starting array for builing container in iterator
    std::array<int, 0> get_atom_ids() const {
      return std::array<int, 0>{};
    }

    //! Access to offsets for access of cluster-related properties
    template <size_t Level, size_t CallerDepth>
    inline size_t get_offset(const ClusterRefBase<Level,
                             CallerDepth> & cluster) const {
      constexpr auto depth{this->cluster_depth<Level>()};
      return cluster.get_cluster_index(depth);
    }

    //! Used for building cluster indices
    template <size_t Level>
    inline size_t get_offset(const std::array<size_t, Level> & counters) const {
      return this->implementation().get_offset_impl(counters);
    }

    inline std::array<size_t, 1> get_counters() const {
      return std::array<size_t, 1>{};
    }

    inline ClusterIndex_t & get_cluster_indices() {
      return this->cluster_indices;
    }

    /**
     * Tuple which contains MaxLevel number of cluster_index lists for reference
     * with increasing depth.  It is filled upon construction of the
     * neighbourhood manager via a
     * std::get<Level>(this->cluster_indices). Higher levels are constructed in
     * adaptors accordingly via the lower level indices and a Level-dependend
     * index is appended to the array.
     */
    ClusterIndex_t cluster_indices;

    // template <size_t Level, size_t Depth>
    // inline int get_cluster(const ClusterRefBase<Level, Depth> & cluster) const {
    //   // all pairs of the following thing
    // }



  private:
  };


  namespace internal {

    template <typename T, size_t Size, int... Indices>
    decltype(auto) append_array_helper(const std::array<T, Size> & arr, T &&  t,
                                       std::integer_sequence<int, Indices...>) {
      return std::array<T, Size+1> {arr[Indices]..., std::forward<T>(t)};
    }

    template <typename T, size_t Size>
    decltype(auto) append_array (const std::array<T, Size> & arr, T &&  t) {
      return append_array_helper(arr, std::forward<T>(t),
                                 std::make_integer_sequence<int, Size>{});
    }

    // template <size_t Level, class AtomRef_t, std::size_t... I>
    // std::array<int, Level>
    // get_indices_from_list(const std::array<AtomRef_t, Level> & atoms,
    //                       std::integer_sequence<int, I...>) {
    //   return std::array<int, Level>{atoms[I].get_index()...};
    // }

    // template <size_t Level, class AtomRef_t> std::array<int, Level>
    // get_indices(const std::array<AtomRef_t, Level> & atoms) {
    //   return get_indices_from_list(atoms, std::make_integer_sequence<int, Level>{});
    // }

    /**
     * Depending on the Level of the ClusterRef (Level=1 or higher), the cluster
     * can be a ghost cluster. If it is a ghost-cluster, i.e. the position of a
     * periodic image will be returned automatically. The decision is made here.
     */
    template <size_t Level, class ClusterRef>
    struct PositionGetter {
      using Vector_ref = typename ClusterRef::Manager_t::Vector_ref;
      static inline Vector_ref get_position(ClusterRef & cluster) {
        return cluster.get_manager().neighbour_position(cluster);
      };
    };

    template <class ClusterRef>
    struct PositionGetter<1, ClusterRef> {
      using Vector_ref = typename ClusterRef::Manager_t::Vector_ref;
      static inline Vector_ref get_position(ClusterRef & cluster) {
        return cluster.get_manager().position(cluster.back());
      };
    };
  }  // internal
  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  class NeighbourhoodManagerBase<ManagerImplementation>::AtomRef
  {
  public:
    using Vector_t = Eigen::Matrix<double, ManagerImplementation::dim(), 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
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

    //! return index
    inline const int & get_index() const {return this->index;}

    //! return position vector
    inline Vector_ref get_position() {
      return this->manager.position(this->index);
    }

    //! return atom type
    inline int get_atom_type() const {
      return this->manager.atom_type(this->index);
    }

  protected:
    Manager_t & manager;
    /**
     * The meaning of `index` is manager-dependent. There are no
     * guaranties regarding contiguity. It is used internally to
     * absolutely address atom-related properties.
     */
    int index;
  private:
  };



  /* ---------------------------------------------------------------------- */
  /**
   *  This is the object we have when iterating over the manager
   */
  template <class ManagerImplementation>
  template <size_t Level>
  class NeighbourhoodManagerBase<ManagerImplementation>::ClusterRef :
    public ClusterRefBase<Level,
                          ManagerImplementation::template cluster_depth<Level>()>
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    using Parent =
      ClusterRefBase<Level,
                     ManagerImplementation::template cluster_depth<Level>()>;
    using AtomRef_t = typename Manager_t::AtomRef;
    using Iterator_t = typename Manager_t::template iterator<Level>;
    using Atoms_t = std::array<AtomRef_t, Level>;
    using iterator = typename Manager_t::template iterator<Level+1>;
    friend iterator;

    using IndexConstArray_t = typename Parent::IndexConstArray;
    using IndexArray_t = typename Parent::IndexArray;

    static_assert(Level <= traits::MaxLevel,
                  "Level > MaxLevel, impossible iterator");

    //! Default constructor
    ClusterRef() = delete;

    //! Constructor from an iterator
    // ClusterRef(Iterator_t & it):
    //   Parent{it.get_atom_indices(), property(it.get_cluster_index())},
    //   it{it} {}

    //! ClusterRef for multiple atoms with const IndexArray
    ClusterRef(Iterator_t & it,
               const std::array<int, Level> & atom_indices,
               const IndexConstArray_t & cluster_indices) :
      Parent{atom_indices, cluster_indices}, it{it} {}

    //! ClusterRef for multiple atoms with non const IndexArray
    ClusterRef(Iterator_t & it,
               const std::array<int, Level> & atom_indices,
               IndexArray_t & cluster_indices) :
      Parent{atom_indices, IndexConstArray_t(cluster_indices.data())}, it{it} {}

    //! ClusterRef for single atom, see `cluster_index`
    ClusterRef(Iterator_t & it,
               const std::array<int, Level> & atom_indices,
               const size_t & cluster_index) :

      Parent{atom_indices, IndexConstArray_t (& cluster_index)}, it{it} {}

    // needed to construct a ClusterRef<1> with the correct offset:


    // Reference to j neighbours of a given atom i // TODO: description
    // cluster.get_indices() -> index of the atom w.r.t order in arrays
    /**
     * This is a ClusterRef of Level=1, constructed from a higher Level.  This
     * function here is self referencing right now. A ClusterRefBase with
     * Level=1 is needed to construct it ?!
     */
    template <bool FirstLevel=(Level==1)>
    ClusterRef(std::enable_if_t<FirstLevel, ClusterRefBase<1,0>> & cluster,
               Manager_t & manager):
      Parent(cluster.get_atom_indices(), cluster.get_cluster_indices()),
      it(manager) {}

    // ClusterRef(Manager_t & manager, size_t access_index);

    //! construct a clusterref from atom_offset: get a clusterref of Level=1 to
    //! iterate over neighbours
    // ClusterRef(const size_t atom_index, Manager_t & manager):
    //   Parent{}

    //! Copy constructor
    ClusterRef(const ClusterRef & other) = delete;

    //! Move constructor
    ClusterRef(ClusterRef && other) = default;

    //! Destructor
    virtual ~ClusterRef() = default;

    //! Copy assignment operator
    ClusterRef& operator=(const ClusterRef & other) = default;

    //! Move assignment operator
    ClusterRef& operator=(ClusterRef && other) = default;


    const std::array<AtomRef_t, Level> & get_atoms() const {
      return this->atoms;
    }

    // TODO: Not sure if this function is needed/used/necessary
    // const std::array<int, Level> & get_atom_ids() const {
    //   return this->atom_indices;
    // }

    const std::array<int, Level> & get_atoms_ids() const {return this->atom_indices;};


    /* There are 2 cases:
     * center (Level== 1)-> position is in the cell
     * neighbour (Level > 1) -> position might have an offset (ghost atom)
     */
    inline Vector_ref get_position() {
      return internal::PositionGetter<Level, ClusterRef>::get_position(*this);
    }

    inline decltype(auto) get_atom_type() {
      auto && id{this->atom_indices.back()};
      return this->get_manager().atom_type(id);
    }

    //! return the index of the atom: Atoms_t is len==1 if center,
    //! len==2 if 1st neighbours,...
    inline int get_atom_index() {
      return this->back(); // TODO: ?? what .back() is it?
    }

    inline Manager_t & get_manager() {return this->it.get_manager();}
    inline const Manager_t & get_manager() const {
      return this->it.get_manager();
    }

    inline iterator begin() {
      std::array<size_t, Level> counters{this->it.get_counters()};
      auto offset = this->get_manager().get_offset(counters);
      return iterator(*this, 0, offset);
    }
    inline iterator end() {
      return iterator(*this, this->size(), std::numeric_limits<size_t>::max());
    }
    inline size_t size() {return this->get_manager().cluster_size(*this);}

    size_t get_index() const {
      return this->it.index;
    }

    inline size_t get_global_index() const {
      return this->get_manager().get_offset(*this);
    }

    const std::array<int, Level> & get_atom_indices() const {
      return this->atom_indices;
    }

  protected:
    inline std::array<size_t, 1> get_counters() const {
      return this->it.get_counters();
    }
    /**
     * `atom_cluster_indices` is an initially contiguous numbering of atoms
     */
    // IndexConstArray_t atom_cluster_indices;
    Iterator_t & it;
  private:
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <size_t Level>
  class NeighbourhoodManagerBase<ManagerImplementation>::iterator
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    friend Manager_t;
    using ClusterRef_t = typename Manager_t::template ClusterRef<Level>;

    friend ClusterRef_t;
    using Container_t =
      std::conditional_t
      <Level == 1,
       Manager_t,
       typename Manager_t::template
       ClusterRef<Level-1>>;
    static_assert(Level > 0, "Level has to be positive");

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

    //! dereference
    //! calculate cluster indices
    inline value_type operator * () {
      auto & cluster_indices_properties =
        std::get<Level-1>(this->get_manager().get_cluster_indices());
      using Ref_t = typename
        std::remove_reference_t<decltype(cluster_indices_properties)>::reference;
      Ref_t cluster_indices =
        cluster_indices_properties[this->get_cluster_index()];
      return ClusterRef_t(*this, this->get_atom_indices(), cluster_indices);
    }

    //! equality
    inline bool operator == (const iterator & other) const {
      return this->index == other.index;
    }

    //! inequality
    inline bool operator != (const iterator & other) const {
      return not (*this == other);
    }

  protected:
    //! constructor with container ref and starting point
    iterator(Container_t & cont, size_t start, size_t offset)
      :container{cont}, index{start}, offset{offset} {}

    //! add atomic indices in current iteration
    std::array<int, Level> get_atom_indices() {
      return internal::append_array
        (container.get_atom_indices(),
         this->get_manager().cluster_neighbour(container, this->index));
    }

    //! returns the current index of the cluster in iteration
    inline size_t get_cluster_index() const {
      return this->index + this->offset;
    }

    inline Manager_t & get_manager() {return this->container.get_manager();}
    inline const Manager_t & get_manager() const {
      return this->container.get_manager();
    }

    inline std::array<size_t, Level> get_counters() {
      std::array<size_t, Level> counters;
      counters[Level-1] = this->index;
      if (Level == 1) {
        return counters;
      } else {
        auto parental_counters = this->container.get_counters();
        for (size_t i{0}; i < Level-1; i++) {
          counters[i] = parental_counters[i];
        }
        return counters;
      }
    }

    Container_t & container;
    size_t index; // + offset from earlier pairs.
    const size_t offset;
  private:
  };

}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_BASE_H */
