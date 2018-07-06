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

#include <Eigen/Dense>

#include <cstddef>
#include <array>
#include <type_traits>
#include <utility>

namespace rascal {



  namespace AdaptorTraits {

    enum class SortedByDistance: bool {yes = true, no = false};
    enum class MinImageConvention: bool {yes = true, no = false};
    enum class NeighbourListType {full, half};
    //----------------------------------------------------------------------------//
    enum class Strict:bool {yes = true, no = false}; // r_cut

    class Type; // type_id
  }  // AdaptorTraits

  //! traits structure to avoid incomplete types in crtp
  template <class Manager>
  struct NeighbourhoodManager_traits
  {};

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
    //! Default constructor
    NeighbourhoodManagerBase() = default;

    //! Copy constructor
    NeighbourhoodManagerBase(const NeighbourhoodManagerBase & other) = delete;

    //! Move constructor
    NeighbourhoodManagerBase(NeighbourhoodManagerBase && other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerBase() = default;

    //! Copy assignment operator
    NeighbourhoodManagerBase & operator=(const NeighbourhoodManagerBase & other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerBase & operator=(NeighbourhoodManagerBase && other)  = default;

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

    inline Iterator_t begin() {return Iterator_t(*this, 0);}
    inline Iterator_t end() {return Iterator_t(*this,
                                               this->implementation().size());}
    inline size_t size() const {return this->implementation().get_size();}

    inline size_t nb_clusters(size_t cluster_size) const {
      return this->implementation().get_nb_clusters(cluster_size);
    }

    inline Vector_ref position(const size_t & atom_index) {
      return this->implementation().get_position(atom_index);
    }
    template <size_t L, size_t D>
    inline Vector_ref neighbour_position(ClusterRefBase<L, D> & cluster) {
      return this->implementation().get_neighbour_position(cluster);
    }

    inline size_t atom_type(const size_t & atom_index) {
      return this->implementation().get_atom_type(atom_index);
    }

  protected:
    template <size_t Level>
    constexpr static size_t cluster_depth(){
      return compute_cluster_depth<Level>(typename traits::DepthByDimension{});
    }

    //! recursion end, not for use
    const std::array<size_t, 0> get_atom_indices() const {
      return std::array<size_t, 0>{};
    }

    template <size_t L>
    inline size_t cluster_size(ClusterRef<L> & cluster) const {
      return this->implementation().get_cluster_size(cluster);
    }

    //! get atom_index of index-th atom in cluster
    template <size_t L>
    inline size_t atom_id(ClusterRef<L> & cluster, size_t index) const {
      return this->implementation().get_atom_id(cluster, index);
    }

    //! get atom_index of the index-th atom in manager
    inline size_t atom_id(NeighbourhoodManagerBase & cluster,
                          size_t & index) const {
      return this->implementation().get_atom_id(cluster, index);
    }

    inline NeighbourhoodManagerBase & get_manager() {return *this;}

    inline ManagerImplementation & implementation() {
      return static_cast<ManagerImplementation&>(*this);
    }
    inline const ManagerImplementation & implementation() const {
      return static_cast<const ManagerImplementation&>(*this);
    }

    std::array<AtomRef, 0> get_atoms() const {return std::array<AtomRef, 0>{};};

    //! Starting array for builing container in iterator
    std::array<size_t, 0> get_atom_ids() const {return std::array<size_t, 0>{};};

    // TODO: get property index - dependent on Depth (e.g., by means of sorting)
    // not necessary to specialize, can be done in _base
    template <size_t Level, size_t CallerDepth>
    inline size_t get_offset(const ClusterRefBase<Level, CallerDepth> & cluster) const {
      constexpr static auto ActiveDepth{
        compute_cluster_depth<Level>(typename traits::DepthByDimension{})};
      static_assert(CallerDepth>=ActiveDepth,
                    "Calling from an inexisting depth");
      return cluster.get_cluster_index(ActiveDepth);
      //      return this->implementation().get_offset_impl(cluster);
    }

    // inline size_t * cluster_indices
    template <size_t Depth>
    inline Eigen::Map<Eigen::Array<size_t, Depth+1, 1>> get_cluster_indices() {
    }

    // template <size_t L>
    // inline int get_cluster(const ClusterRefBase<L> & cluster) const {
    //   // all pairs of the following thing
    // }

  private:
  };


  namespace internal {
    
    // OLD VERSION ----------
    template <typename T, size_t Size, size_t... Indices>
    decltype(auto) append_array_helper(const std::array<T, Size> & arr, T &&  t,
                                        std::index_sequence<Indices...>) {
      return std::array<T, Size+1> {arr[Indices]..., std::forward<T>(t)};
    }
    
    template <typename T, size_t Size>
    decltype(auto) append_array (const std::array<T, Size> & arr, T &&  t) {
      return append_array_helper(arr, std::forward<T>(t),
                                 std::make_index_sequence<Size>{});
    }

    // CHANGED VERSION ----
    // template <typename T, size_t Size, size_t... Indices>
    // decltype(auto) append_array_helper(std::array<T, Size> && arr, T &&  t,
    //                                    std::index_sequence<Indices...>) {
    //   return std::array<T, Size+1> {(arr[Indices])..., std::forward<T>(t)};
    // }

    // // TODO: not corresponding to
    // //append_array(std::remove_reference<const std::array<long unsigned int, 0>&>::type, size_t)
    // //std::conditional_t<Level == 1, Manager_t, typename Manager_t::template ClusterRef<Level-1>>;
    // template <typename T, size_t Size>
    // decltype(auto) append_array(std::array<T, Size> && arr, T &&  t) {
    //   return append_array_helper(arr, std::forward<T>(t),
    //                              std::make_index_sequence<Size>{});
    // }

    template<size_t Level, class AtomRef_t, std::size_t... I>
    std::array<size_t, Level>
    get_indices_from_list(const std::array<AtomRef_t, Level> & atoms,
                          std::index_sequence<I...>) {
      return std::array<size_t, Level>{atoms[I].get_index()...};
    }

    template<size_t Level, class AtomRef_t> std::array<size_t, Level>
    get_indices(const std::array<AtomRef_t, Level> & atoms) {
      return get_indices_from_list(atoms, std::make_index_sequence<Level>{});
    }

    template<size_t Level, class ClusterRef>
    struct PositionGetter {
      using Vector_ref = typename ClusterRef::Manager_t::Vector_ref;
      static inline Vector_ref  get_position(ClusterRef & cluster) {
        return cluster.get_manager().neighbour_position(cluster);
      };
    };

    template<class ClusterRef>
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
    //AtomRef(Manager_t & manager, int id): manager{manager}, index{id}{}
    AtomRef(Manager_t & manager, const size_t & id): manager{manager}, index{id} {}
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
    inline const size_t & get_index() const {return this->index;}

    //! return position vector
    inline Vector_ref get_position() {return this->manager.position(this->index);}

    //! return atom type
    inline size_t get_atom_type() const {return this->manager.atom_type(this->index);}

  protected:
    Manager_t & manager;
    const size_t & index; //!< corresponds to the 0-depth, 1st-level cluster_id (i.e. memory-offset)
  private:
  };



  /* ---------------------------------------------------------------------- */
  /**
    This is the object we have when iterating over the manager
  */
  template <class ManagerImplementation>
  template <size_t Level>
  class NeighbourhoodManagerBase<ManagerImplementation>::ClusterRef :
    public ClusterRefBase<Level, ManagerImplementation::template cluster_depth<Level>()>
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    using Parent = ClusterRefBase<Level, ManagerImplementation::template cluster_depth<Level>()>;
    using AtomRef_t = typename Manager_t::AtomRef;
    using Iterator_t = typename Manager_t::template iterator<Level>;
    using Atoms_t = std::array<AtomRef_t, Level>;
    using iterator = typename Manager_t::template iterator<Level+1>;
    friend iterator;

    static_assert(Level <= traits::MaxLevel,
                  "Level > MaxLevel, impossible iterator");

    //! Default constructor
    ClusterRef() = delete;

    //! Constructor from an iterator
    ClusterRef(Iterator_t & it):
      Parent{it.get_atom_indices(), it.get_cluster_indices()},
      //it{it}{}
        // Intelligenzija: an array of cluster_indeces
        // std::array<size_t, Depth>
      it{it}
    {}


    template<size_t Depth>
    ClusterRef(std::enable_if<Level==1, ClusterRefBase<1, Depth>> & cluster,
               Manager_t& manager):
      Parent{cluster.get_indices(), cluster.get_cluster_indices()},
      it{manager}{}

    //! Copy constructor
    ClusterRef(const ClusterRef & other) = default;

    //! Move constructor
    ClusterRef(ClusterRef && other) = default;

    //! Destructor
    virtual ~ClusterRef() = default;

    //! Copy assignment operator
    ClusterRef& operator=(const ClusterRef & other) = default;

    //! Move assignment operator
    ClusterRef& operator=(ClusterRef && other) = default;


    const std::array<AtomRef_t, Level> & get_atoms() const {return this->atoms;};
    std::array<AtomRef_t, Level> & get_atoms() {return this->atoms;};

    const std::array<size_t, Level> & get_atom_ids() const {return this->atom_indices;};
    std::array<size_t, Level> & get_atoms_ids() {return this->atom_indices;};


    /* There are 2 cases:
        center (Level== 1)-> position is in the cell
        neighbour (Level > 1)   -> position might have an offset associated
     */

    inline Vector_ref get_position() {
      return internal::PositionGetter<Level, ClusterRef>::get_position(*this);
    }

    inline decltype(auto) get_atom_type() const {
      return this->atoms.back().get_atom_type();
    }

    //! return the index of the atom: Atoms_t is len==1 if center,
    //! len==2 if 1st neighbours,...
    inline size_t get_atom_index() {
      return this->back();
    }

    inline Manager_t & get_manager() {return this->it.get_manager();}
    inline const Manager_t & get_manager() const {return this->it.get_manager();}

    inline iterator begin() {return iterator(*this, 0);}
    inline iterator end() {return iterator(*this, this->size());}
    inline size_t size() {return this->get_manager().cluster_size(*this);}
    size_t get_index() const {
      return this->it.index;
    }

    inline size_t get_global_index() const {
      return this->get_manager().get_offset(*this);
    }

  protected:
    //Atoms_t atoms; // TODO no atom refs any more, they are in cluster_ref_base now
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
    inline value_type operator * () {
      return ClusterRef_t(*this);
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
    iterator(Container_t & cont, size_t start)
      :container{cont}, index{start} {}

    std::array<size_t, Level> get_atom_indices() {
      static_assert(std::is_same<size_t, long unsigned int>::value,
                    "mhm");
      return internal::append_array
        (container.get_atom_indices(),
         this->get_manager().atom_id(container, this->index));
    }

    // TODO: something like this for array size
    // this->get_manager().compute_cluster_depth<Level>(typename traits::DepthByDimenions{})+1
    Eigen::Map<
      Eigen::Array<size_t,
                   compute_cluster_depth<Level>(typename traits::DepthByDimension{})+1,
                // ActiveDepth+1,
               // should be doing something like this: DEPTH+1,
                   1>> get_cluster_indices() {
      // TODO: global index from offset?
      // needs its own global index
      // this array need to be filled with this->index upon iteration
      return Eigen::Map<Eigen::Array<size_t, compute_cluster_depth<Level>(typename traits::DepthByDimension{})+1, 1>>(nullptr);
    }

    // TODO: This is obsolete, because the atoms are not any more a part of the cluster!
    // std::array<AtomRef_t, Level> get_container_atoms() {
    //   return internal::append_array
    //     (std::move(container.get_atoms()),  // TODO:  should be container.get_atom_indeces?
    //      AtomRef_t(this->get_manager(),
    //                this->get_manager().atom_id(container, this->index)));
    // }

    inline Manager_t & get_manager() {return this->container.get_manager();}
    inline const Manager_t & get_manager() const {
      return this->container.get_manager();
    }

    Container_t & container;
    size_t index;
    std::array<size_t, traits::MaxLevel-Level> starting_points{}; // TODO: initialized by iteration?
  private:
  };


}  // rascal


#endif /* NEIGHBOURHOOD_MANAGER_BASE_H */
