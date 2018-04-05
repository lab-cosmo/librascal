/**
 * file   neighbourhood_manager_base.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Interface for neighbourhood managers
 *
 * @section LICENSE
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
     * iterator over the atoms in the manager. Iterators like these
     * can be used as indices for random access in atom-related
     * fields.
     */
    class iterator;

    /**
     * return type for iterators: a light-weight atom reference,
     * giving access to an atom's position and force
     */
    class AtomRef;

    inline iterator begin() {return iterator(*this, 0);}
    inline iterator end() {return iterator(*this, this->implementation().size());}
    inline size_t size() const {return this->implementation().get_size();}

    inline Vector_ref get_x(const AtomRef atom);

    inline Vector_ref get_f(const AtomRef atom);
  protected:
    inline ManagerImplementation& implementation() {
      return static_cast<ManagerImplementation&>(*this);
    }
    inline const ManagerImplementation& implementation() const {
      return static_cast<const ManagerImplementation&>(*this);
    }
  private:
  };

  template <class ManagerImplementation>
  class NeighbourhoodManagerBase<ManagerImplementation>::AtomRef
  {
  public:
    using Vector_t = Eigen::Matrix<double, ManagerImplementation::dim(), 1>;
    using Vector_ref = Eigen::Map<Vector_t>;
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    using iterator_t = typename Manager_t::iterator;

    //! Default constructor
    AtomRef() = delete;

    //! constructor from iterator
    AtomRef(iterator_t & it): it{it}{}

    //! Copy constructor
    AtomRef(const AtomRef &other) = default;

    //! Move constructor
    AtomRef(AtomRef &&other) = default;

    //! Destructor
    virtual ~AtomRef() = default;

    //! Copy assignment operator
    AtomRef& operator=(const AtomRef &other) = default;

    //! Move assignment operator
    AtomRef& operator=(AtomRef &&other) = default;

    //! return index
    inline int get_index() const {return this->it.index;}

    //! return position vector
    inline Vector_ref get_x() {return this->it.manager.get_x(*this);}
    //! return force vector
    inline Vector_ref get_f() {return this->it.manager.get_f(*this);}

  protected:
    iterator_t & it;
  private:
  };


  template <class ManagerImplementation>
  class NeighbourhoodManagerBase<ManagerImplementation>::iterator
  {
  public:
    using Manager_t = NeighbourhoodManagerBase<ManagerImplementation>;
    friend Manager_t;
    using AtomRef_t = typename Manager_t::AtomRef;
    friend AtomRef_t;

    using value_type = AtomRef_t;
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
      return AtomRef(*this);
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
    //! constructor with manager and starting point
    iterator(Manager_t & manager, int start)
      :manager{manager}, index{start} {}

    Manager_t & manager;
    int index;
  private:
  };


  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  auto
  NeighbourhoodManagerBase<ManagerImplementation>::
  get_x(const AtomRef atom) -> Vector_ref {
    return this->implementation().get_x(std::move(atom));
  }

  /* ---------------------------------------------------------------------- */
  template <class ManagerImplementation>
  auto
  NeighbourhoodManagerBase<ManagerImplementation>::
  get_f(const AtomRef atom) -> Vector_ref {
    return this->implementation().get_f(std::move(atom));
  }

}  // proteus


#endif /* NEIGHBOURHOOD_MANAGER_BASE_H */
