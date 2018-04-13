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
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
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

#include <Eigen/Dense>

#include <cstddef>
#include <array>
#include <type_traits>


namespace proteus {

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
    using Vector_block = Eigen::Block<Eigen::MatrixXd, -1, 1, true>;
    //! Default constructor
    NeighbourhoodManagerBase() = default;

    //! Copy constructor
    NeighbourhoodManagerBase(const NeighbourhoodManagerBase &other) = delete;

    //! Move constructor
    NeighbourhoodManagerBase(NeighbourhoodManagerBase &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerBase() = default;

    //! Copy assignment operator
    NeighbourhoodManagerBase& operator=(const NeighbourhoodManagerBase &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerBase& operator=(NeighbourhoodManagerBase &&other)  = default;

    /**
     * iterator over the atoms, pairs, triplets, etc in the
     * manager. Iterators like these can be used as indices for random
     * access in atom-, pair, ... -related fields.
     */
    template <int Level, int MaxLevel>
    class iterator;
    using Iterator_t = iterator<1, traits::MaxLevel>;
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
    template <int Level, int MaxLevel>
    class ClusterRef;

    inline Iterator_t begin() {return Iterator_t(*this, 0);}
    inline Iterator_t end() {return Iterator_t(*this,
                                               this->implementation().size());}
    inline size_t size() const {return this->implementation().get_size();}

    inline size_t nb_clusters(int cluster_size) const {
      return this->implementation().get_nb_clusters(cluster_size);
    }
    /*
    inline Vector_block get_position(const AtomRef& atom) {
      return this->implementation().get_position(atom);
    }
    */
   inline Vector_block get_position(const AtomRef& atom) {
      return this->implementation().get_position(atom);
    }
    inline Vector_block get_f(const AtomRef& atom) {
      return this->implementation().get_f(atom);
    }

  protected:
    template <int L, int ML>
    inline size_t cluster_size(ClusterRef<L, ML> & cluster) const {
      return this->implementation().get_cluster_size(cluster);
    }
    template <int L, int ML>
    inline size_t atom_id(ClusterRef<L, ML> & cluster, int index) const {
      return this->implementation().get_atom_id(cluster, index);
    }

    inline size_t atom_id(NeighbourhoodManagerBase & cluster, int index) const {
      return this->implementation().get_atom_id(cluster, index);
    }

    inline NeighbourhoodManagerBase & get_manager() {return *this;}

    inline ManagerImplementation& implementation() {
      return static_cast<ManagerImplementation&>(*this);
    }
    inline const ManagerImplementation& implementation() const {
      return static_cast<const ManagerImplementation&>(*this);
    }

    std::array<AtomRef, 0> get_atoms() const {return std::array<AtomRef, 0>{};};

    template <int L, int ML>
    inline int get_offset(const ClusterRef<L, ML> & cluster) const {
      return this->implementation().get_offset_impl(cluster);
    }

  private:
  };


  namespace internal {

    template <typename T, size_t Size, size_t... Indices>
    decltype(auto) append_array_helper(std::array<T, Size>&& arr, T &&  t,
                                        std::index_sequence<Indices...>) {
      return std::array<T, Size+1> {std::move(arr[Indices])..., std::forward<T>(t)};
    }
    template <typename T, size_t Size>
    decltype(auto) append_array (std::array<T, Size>&& arr, T &&  t) {
      return append_array_helper(std::move(arr), std::forward<T>(t),
                                 std::make_index_sequence<Size>{});
    }

  }  // internal
  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  class NeighbourhoodManagerBase<ManagerImplementation>::AtomRef
  {
  public:
    using Vector_t = Eigen::Matrix<double, ManagerImplementation::dim(), 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using Vector_block = Eigen::Block<Eigen::MatrixXd, -1, 1, true>;
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;

    //! Default constructor
    AtomRef() = delete;

    //! constructor from iterator
    AtomRef(Manager_t & manager, int id): manager{manager}, index{id}{}

    //! Copy constructor
    AtomRef(const AtomRef &other) = default;

    //! Move constructor
    AtomRef(AtomRef &&other) = default;

    //! Destructor
    virtual ~AtomRef() = default;

    //! Copy assignment operator
    AtomRef& operator=(const AtomRef &other) = delete;

    //! Move assignment operator
    AtomRef& operator=(AtomRef &&other) = default;

    //! return index
    inline int get_index() const {return this->index;}

    //! return position vector
    // inline Vector_block get_position() {return this->manager.get_position(*this);}
    //! return position vector
    inline Vector_block get_position() {return this->manager.get_position(*this);}
    //! return force vector
    inline Vector_block get_f() {return this->manager.get_f(*this);}

  protected:
    Manager_t & manager;
    int index;
  private:
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <int Level, int MaxLevel>
  class NeighbourhoodManagerBase<ManagerImplementation>::ClusterRef
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    using AtomRef_t = typename Manager_t::AtomRef;
    using Iterator_t = typename Manager_t::template iterator<Level, MaxLevel>;
    using Atoms_t = std::array<AtomRef_t, Level>;

    using iterator = typename Manager_t::template iterator<Level + 1, MaxLevel>;
    friend iterator;

    //! Default constructor
    ClusterRef() = delete;

    //! constructor from an iterator
    ClusterRef(Iterator_t & it): atoms{it.get_container_atoms()}, it{it}{}

    //! Copy constructor
    ClusterRef(const ClusterRef &other) = default;

    //! Move constructor
    ClusterRef(ClusterRef &&other) = default;

    //! Destructor
    virtual ~ClusterRef() = default;

    //! Copy assignment operator
    ClusterRef& operator=(const ClusterRef &other) = default;

    //! Move assignment operator
    ClusterRef& operator=(ClusterRef &&other) = default;


    const std::array<AtomRef_t, Level>& get_atoms() const {return this->atoms;};
    std::array<AtomRef_t, Level>& get_atoms() {return this->atoms;};

    /**
     * convenience functions, because in loops, we frequently like to
     * use clusters as proxies to their last atom
     */
    inline decltype(auto) get_position() {return this->atoms.back().get_position();}
    inline decltype(auto) get_f() {return this->atoms.back().get_f();}

    inline Manager_t & get_manager() {return this->it.get_manager();}
    inline const Manager_t & get_manager() const {return this->it.get_manager();}

    inline iterator begin() {return iterator(*this, 0);}
    inline iterator end() {return iterator(*this, this->size());}
    inline size_t size() {return this->get_manager().cluster_size(*this);}
    inline int get_index() const {
      return this->it.index;
    }
    inline int get_global_index() const {
      return this->get_manager().get_offset(*this);
    }
  protected:
    Atoms_t atoms;
    Iterator_t & it;
  private:
  };

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  template <int Level, int MaxLevel>
  class NeighbourhoodManagerBase<ManagerImplementation>::iterator
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    friend Manager_t;
    using ClusterRef_t = typename Manager_t::template ClusterRef<Level, MaxLevel>;
    friend ClusterRef_t;
    using Container_t =
      std::conditional_t
      <Level == 1,
       Manager_t,
       typename Manager_t::template ClusterRef<Level-1, MaxLevel>>;
    static_assert(Level > 0, "Level has to be positive");

    using AtomRef_t = typename Manager_t::AtomRef;

    using value_type = ClusterRef_t;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    using reference = value_type;

    //! Default constructor
    iterator() = delete;

    //! Copy constructor
    iterator(const iterator &other) = default;

    //! Move constructor
    iterator(iterator &&other) = default;

    //! Destructor
    virtual ~iterator() = default;

    //! Copy assignment operator
    iterator& operator=(const iterator &other) = default;

    //! Move assignment operator
    iterator& operator=(iterator &&other) = default;

    //! pre-increment
    inline iterator & operator ++ () {
      ++this->index;
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
    iterator(Container_t & cont, int start)
      :container{cont}, index{start} {}

    std::array<AtomRef_t, Level> get_container_atoms() {
      return internal::append_array
        (std::move(container.get_atoms()),
         AtomRef_t(this->get_manager(),
                   this->get_manager().atom_id(container, this->index)));
    }

    inline Manager_t & get_manager() {return this->container.get_manager();}
    inline const Manager_t & get_manager() const {
      return this->container.get_manager();
    }

    Container_t & container;
    int index;
  private:
  };


}  // proteus


#endif /* NEIGHBOURHOOD_MANAGER_BASE_H */
